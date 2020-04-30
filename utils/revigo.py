#!/usr/bin/python

"""
- Submit example data to REVIGO server (http://revigo.irb.hr/)
- Download and run R script for creating the treemap
- Download and run R script for creating the scatterplot

Creates files:
    treemap.R, treemap.Rout, revigo_treemap.pdf
    scatter.R, scatter.Rout, revigo_scatter.pdf
"""

import os, sys
import urllib
import mechanize
import shutil
import pandas as pd
import numpy as np

args = sys.argv


def sumbit_to_revigo(base_folder="/home/hag007/Desktop/aggregate_gwas_report/oob",file_format="emp_diff_modules_{}_{}_passed_oob.tsv", dataset="Breast_Cancer.G50", algorithm="my_netbox_td"):
    fl=os.path.join(base_folder, file_format.format(dataset, algorithm))
    os.chdir(base_folder)

    url = "http://revigo.irb.hr/"

    # RobustFactory because REVIGO forms not well-formatted
    br = mechanize.Browser()

    # For actual data, use open('mydata.txt').read()
    # br.open(os.path.join(url, 'examples', 'example1.txt'))
    # txt = br.response().read()
    mtx = pd.read_csv(fl, error_bad_lines=False, sep='\t').dropna(axis=0)[['GO id','hg_pval_max']]
    mtx["hg_pval"]=mtx['hg_pval_max'].apply(lambda a: 10**(-a))
    np.savetxt("foo.csv", mtx[['GO id', 'hg_pval']], delimiter="\t", fmt='%s')
    txt=open("foo.csv").read()
    # txt="GO:0009268 1e-14\nGO:0010447 1e-14"
    # directory = args[1].split(".")[0]
    # directory = "../revigo/" + args[2] + "/" + directory
    # # shutil.rmtree(directory)
    # os.makedirs(directory)

    # Encode and request
    data = {'inputGoList': txt}
    br.open(url, data=urllib.urlencode(data))

    # Choose HomoSapiens in the drop down menu
    br.select_form(nr=0)
    control = br.form.find_control("goSizes")
    # loop through items to find the match
    for item in control.items:
        if item.name == "9606":
            # it matches, so select it
            item.selected = True

    control = br.form.find_control("cutoff")
    # loop through items to find the match
    for item in control.items:
        if item.name == "0.50":
            # it matches, so select it
            item.selected = True

    control = br.form.find_control("measure")
    # loop through items to find the match
    for item in control.items:
        if item.name == "RESNIK":
            # it matches, so select it
            item.selected = True

    # Submit form
    br.select_form(name="submitToRevigo")
    response = br.submit()

    # Exact string match on the url for getting the R treemap script
    br.follow_link(url="export.jsp?table=1")
    with open('REVIGO.csv', 'w') as f:
        f.write(br.response().read())

    br.back()

    br.follow_link(url="toR_treemap.jsp?table=1")
    with open('treemap.R', 'w') as f:
        f.write(br.response().read())

    # Create PDFs
    # os.chdir(args[2])
    os.system('Rscript treemap.R')

if __name__=="__main__":
    sumbit_to_revigo()