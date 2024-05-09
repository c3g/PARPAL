import os
import sys
from shiny import ui, render, reactive, App
import shinyswatch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from matplotlib import colors
import asyncio

## root of data  
DATAROOT=os.getenv("PARPAL_DATA_LOCATION", "/data")

## Load data and compute static values
meta=pd.read_table(f"{DATAROOT}/metadata.tsv")## CHANGE HERE ##

## add DATAROOT env var to gfp path.   
meta["gfp path"]= f"{DATAROOT}/" + meta['gfp path']

## Select genes of interes to choose from 
genesofi=sorted(meta["pairs display"].unique().tolist())

## Select my picture colors
mycmap=colors.LinearSegmentedColormap.from_list("custom", colors=["#000000",'#83f52c'], N=50)


################################################################################
##################### ASSUMPTIONS FOR THIS TO WORK #############################

# 0 - Metadata must be in /data folder.  
# 1 - There must always be at 3 replicates, no more no less 
# 2 - There are only 2 Tpairs for each gene pair
# 3 - Only one gene-pair at the time
# 4 - Assuming "label" is always the same in the same format.  

##################### ASSUMPTIONS FOR THIS TO WORK #############################
################################################################################


# app_ui #######################################################################
app_ui = ui.page_sidebar(
    ## Beginnning of Sidebar 
    ui.sidebar(
        ui.input_selectize(id = "genepairs", label = "Search for Paralogs:"
          , choices= genesofi
          , multiple = False
          , selected = None 
          , remove_button = True
        ).add_style("font-family:Helvetica Neue;"),
        # ui.h1("Test"),
        ui.input_action_button("submit", "Submit", class_="btn-success").add_style("background:#7393B3; border:white; font-family:Helvetica Neue;"), 
        ui.output_ui("compute"),
        ui.output_text_verbatim("C3G", placeholder = False).add_style("text-align:center; background:#F8F8F8; border:white; font-size:9pt; font-family:Helvetica Neue;"),
        ui.output_image("C3Gimage", inline = True).add_style("text-align:center; padding:10px;"),
        shinyswatch.theme.litera(),
        width=400,
        padding=40, ## Padding btw things inside
    ),
    ## Beginning of the main page
    # ui.h1("Test"),
    ui.output_image("EKimage", inline = True).add_style("text-align:center; padding-right:100px; padding-left:100px;"),
    ui.h3("Introduction").add_style("text-align:center; padding:40px; font-family:Helvetica Neue; font-size:22pt;"),
    ui.div(), 
    ui.output_text_verbatim("message", placeholder = True).add_style("text-align:center; padding:10px; background:white; font-size:14pt; border:white; font-family:Helvetica Neue;"),
    ui.page_navbar(
        ui.nav_panel("Scores", 
            ui.output_table("paralogabund").add_style("text-align:center; padding-top:50px; padding-right:300px; padding-left:300px; font-size:11pt; font-family:Helvetica Neue;"),
            ui.output_table("paralogredis").add_style("text-align:center; padding-top:30px; padding-bottom:50px; padding-right:300px; padding-left:300px; font-size:10pt; font-family:Helvetica Neue;"),
        ),
        ui.nav_panel("Replicate 1",
            ui.output_plot("rep1", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:50px;"),
        ),
        ui.nav_panel("Replicate 2", 
            ui.output_plot("rep2", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:50px;"),
        ),
        ui.nav_panel("Replicate 3", 
            ui.output_plot("rep3", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:50px;"),
        ),
        id="tab",
    ),
    ui.h3("References").add_style("text-align:center; padding:40px; font-family:Helvetica Neue;font-size:22pt;"),
    ui.output_text_verbatim("acknowledge").add_style("text-align:center; padding:10px; background:white; font-size:14pt; border:white; font-family:Helvetica Neue;"),
    # shinyswatch.theme.lux(),## Theme
    # shinyswatch.theme.minty(),
    shinyswatch.theme.litera(),
    title="PARPAL",## Web title
)

# Server ######################################################################
def server(input, output, session):
  
    ## Add Image for side bar
    @render.image
    def C3Gimage():
        img: ImgData = {"src": f"{DATAROOT}/images/c3g.jpg", "width": "50%"}
        return img
    
    ## Set code block for sidebar
    @render.text
    def C3G():
        return f"\n\n\n\n\n\n\n\n\n\n\n\n" + \
                f"Website developed by: \n" + \
                f"Gerardo Zapata, Rohan Dandage, \n" + \
                f"Vanessa Pereira and Elena Kuzmin \n\n" + \
                f"In Collaboraton with: \n" + \
                f"Canadian Centre for Computational Genomics (C3G) \n" + \
                f"https://computationalgenomics.ca/team_profiles/#GZ"

    ## Make little loading feature - form submit in sidebar
    @output
    @render.ui
    @reactive.event(input.submit)
    async def compute():
        with ui.Progress(min=1, max=95) as p:
            p.set(message="Calculation in progress", detail="This may take a while...")

            for i in range(1, 35):
                p.set(i, message="Computing")
                await asyncio.sleep(0.05)

        return

    ## Add Image for main page
    @render.image
    def EKimage():
        img: ImgData = {"src": f"{DATAROOT}/images/20240417_PARPAL_logo.png", "width": "80%"}
        return img
    
    ## remake metadata table --> for gene pair
    @reactive.calc
    def meta2():
        meta2=meta[meta["pairs display"] == input.genepairs()].sort_values(by='label', ascending=False).reset_index(drop = True)
        return meta2

    ## make variable for the two genes 
    @reactive.calc
    def gene1():
        meta_of_step=meta2()
        genes_of_step=meta_of_step["pairs"][1].split("-")
        gene1=genes_of_step[0]
        return gene1

    @reactive.calc
    def gene2():
        meta_of_step=meta2()
        genes_of_step=meta_of_step["pairs"][1].split("-")
        gene2=genes_of_step[1]
        return gene2
            
    ## Set code block for text example
    @render.text
    def message():
        return f"PARPAL is a web database for single-cell imaging \n" + \
                f"of protein dynamics of paralogs revealing mechanisms of gene retention \n\n" + \
                f"Database statistics:  \n" + \
                f"Proteins screened = 164 \n" + \
                f"Paralog pairs screened = 82 \n" + \
                f"Total micrographs = ~3.5K \n" + \
                f"Total cells = ~460K \n\n" + \
                f"For details on the PARPAL project, data and website please contact Elena Kuzmin: \n" + \
                f"Email: elena.kuzmin@concordia.ca \n" + \
                f"Website: https://kuzmin-lab.github.io/ \n\n"

    ## abundance table
    @render.table
    @reactive.event(input.submit)
    def paralogabund():
        meta_of_step=meta2()
        
        ## Set table options for 3 decimal points 
        pd.options.display.float_format = "{:,.3f}".format 
        
        abund=pd.read_csv(f"{DATAROOT}/scores/by_label_abundance.tsv", sep='\t', float_precision='round_trip')
        abund=abund[abund['label'] == meta_of_step["label"][0]]
        abund=abund[['Paralog pair','protein abundance mean']]
        
        return abund
                
    ## redistribution table 
    @render.table
    @reactive.event(input.submit)
    def paralogredis():
        gene_of_step1=gene1()
        gene_of_step2=gene2()
        
        ## Set table options for 3 decimal points 
        pd.options.display.float_format = "{:,.3f}".format 
        
        redis=pd.read_csv(f"{DATAROOT}/scores/by_gene_redistribution.tsv", sep='\t', float_precision='round_trip')
        redis=redis[(redis['gene symbol'] == gene_of_step1) | (redis['gene symbol'] == gene_of_step2)]
        
        return redis

    ## Set variable for max intensity value per replicate 
    @reactive.calc
    @reactive.event(input.submit)
    def maxs1():
        maxs=[0]
        return maxs

    ## Set variable for max intensity value per replicate 
    @reactive.calc
    @reactive.event(input.submit)
    def maxs2():
        maxs=[0]
        return maxs

    ## Set variable for max intensity value per replicate 
    @reactive.calc
    @reactive.event(input.submit)
    def maxs3():
        maxs=[0]
        return maxs
          
    
    ## Set text for replicate 1 ## CHANGE HERE ##
    @render.text
    @reactive.event(input.submit)
    def rep1t(): ## CHANGE HERE ##
        return f"This is the Beginning of Replicate 1"## CHANGE HERE ##
      
    ## Set plot for replicate 1 ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def rep1():## CHANGE HERE ##
        meta_of_step=meta2()
        
        if len(meta_of_step[meta_of_step["replicate"] == "replicate1"]) > 0:## CHANGE HERE ##
            
            ## set meta for replicate 1
            meta_of_step=meta_of_step[meta_of_step["replicate"] == "replicate1"].sort_values(by='label', ascending=False).reset_index(drop = True)## CHANGE HERE ## 
            
            ## set gfp path for range
            gfp_path=meta_of_step["gfp path"]

            ## Get range of both labels
            ima_range=range(meta_of_step.index.min()
                          , meta_of_step.index.max() + 1
                          )
                          
            ## Get maxs limits
            maxlist=maxs1()
            for i in ima_range:
                im1=imread(fname = gfp_path[i])
                maxi=im1.max()
                maxlist.append(maxi)
                maxlistT=maxlist[1:]
            
            ## T n figures
            n_figs=max(ima_range)+1

            ## Set number of rows based on # of labels 
            n_rows=len(meta_of_step["label"].unique())
            
            if n_rows > 2 :
            
                ## Set fig --> for entire set of Axs
                fig = plt.figure(figsize=(10,10), layout = 'tight')
                
                ## Set number of rows based on # of labels 
                subplots = fig.subfigures(n_rows, 1)
                
                #### Set Loop per number of labels
                for j in range(0, n_rows):
                  
                    ## n figures per row
                    nfigsr=round(n_figs/n_rows)
        
                    ## range of nfigsr  --> n  *** when there are exactly 3 rows
                    if len(meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]]) > 4 and j < 2 :
                        range3=range(0,nfigsr*2+1)
                    elif len(meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]]) < 5 and j == 0 :
                        range3=range(0, nfigsr+1)
                    else : 
                        range3=range(nfigsr,n_figs+1)        
                    # range3=range(nfigsr*j,nfigsr*(j+1)+1)
                    
                    ## range *** when there are exactly 4 rows
                    range4=range(0, nfigsr*2+1) if j < 2 else range(nfigsr*2, n_figs+1)
                    
                    ## set maxlistTT
                    maxlistTT=np.median(maxlistT)*0.3 if n_rows < 3 else min(maxlistT[min(range3):max(range3)]) if n_rows == 3 else min(maxlistT[min(range4):max(range4)])

                    ## If maxlistTT is less than then double
                    maxlistTT=maxlistTT*0.5 if maxlistTT < 300 and n_rows > 2 else maxlistTT
                    
                    ## maxlistTT only for HTA1-HTA2 -- UGLY
                    maxlistTT=maxlistTT*0.55 if maxlistTT < 200 and n_rows > 2 and meta_of_step["pairs display"][0] == "HTA1-HTA2" else maxlistTT*0.2 if maxlistTT > 1000 and n_rows > 2 and meta_of_step["pairs display"][0] == "HTA1-HTA2" else maxlistTT
                
                    ## number of figs per "label" --> for the first unique label 
                    n_cols=len(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]])
                    
                    locals()["ax" + str(j)] = subplots[j].subplots(1, n_cols, sharex = True, sharey=True, squeeze = False) if n_cols == 1 else subplots[j].subplots(1, n_cols, sharex = True, sharey=True)
                    locals()["ax" + str(j)] = locals()["ax" + str(j)].ravel() 
                    
                    for i in range(0, n_cols):
                        test=imread(fname = meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['gfp path'].reset_index(drop = True)[i])
                        locals()["ax" + str(j)][i].imshow(X=test
                                      , cmap = mycmap
                                      , vmax = maxlistTT
                                      )
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['label display'].reset_index(drop = True)[i], size=8) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].set_xticks([])
                        locals()["ax" + str(j)][i].set_yticks([])
                
                fig.suptitle('   ', size=50)
                return fig
            
            else:
                j=0
                
                ## number of figs per "label" --> for the first unique label 
                n_cols=len(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]])##CHANGE HERE##
                
                fig, axs = plt.subplots(1, n_cols, sharex = True, sharey=True, squeeze = False, figsize=(5, 20))
                axs = axs.ravel()
                
                ## Loop per axs
                for i in range(0, n_cols):
                    test=imread(fname = meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['gfp path'].reset_index(drop = True)[i])
                    axs[i].imshow(X=test, cmap = mycmap
                                  , vmax = np.median(maxlistT)*0.3
                                  )
                    axs[i].set_title(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['label display'].reset_index(drop = True)[i], size=8)
                    axs[i].set_xticks([])
                    axs[i].set_yticks([])
                        
                fig.suptitle('   ', size=50)
                return fig

        else:
          
            ## Notification "figure" for no image replicate
            fig, ax = plt.subplots()
            ax.axis([0, 1, 0, 1])
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.spines['bottom'].set_color('white')
            ax.spines['top'].set_color('white') 
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.set_title('No Replicate 1', size=15)## CHANGE HERE ##

            return fig

    ## Set plot for replicate 2 ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def rep2():## CHANGE HERE ##
        meta_of_step=meta2()
        
        if len(meta_of_step[meta_of_step["replicate"] == "replicate2"]) > 0:## CHANGE HERE ##
            
            ## set meta for replicate 1
            meta_of_step=meta_of_step[meta_of_step["replicate"] == "replicate2"].sort_values(by='label', ascending=False).reset_index(drop = True)## CHANGE HERE ## 
            
            ## set gfp path for range
            gfp_path=meta_of_step["gfp path"]

            ## Get range of both labels
            ima_range=range(meta_of_step.index.min()
                          , meta_of_step.index.max() + 1
                          )
                          
            ## Get maxs limits
            maxlist=maxs2()
            for i in ima_range:
                im1=imread(fname = gfp_path[i])
                maxi=im1.max()
                maxlist.append(maxi)
                maxlistT=maxlist[1:]
            
            ## T n figures
            n_figs=max(ima_range)+1
            
            ## Set number of rows based on # of labels 
            n_rows=len(meta_of_step["label"].unique())
            
            if n_rows > 2 :
            
                ## Set fig --> for entire set of Axs
                fig = plt.figure(figsize=(10,10), layout = 'tight')
                
                ## Set number of rows based on # of labels 
                subplots = fig.subfigures(n_rows, 1)
                
                #### Set Loop per number of labels
                for j in range(0, n_rows):
                  
                    ## n figures per row
                    nfigsr=round(n_figs/n_rows)
        
                    ## range of nfigsr  --> n  *** when there are exactly 3 rows
                    if len(meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]]) > 4 and j < 2 :
                        range3=range(0,nfigsr*2+1)
                    elif len(meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]]) < 5 and j == 0 :
                        range3=range(0, nfigsr+1)
                    else : 
                        range3=range(nfigsr,n_figs+1)        
                    # range3=range(nfigsr*j,nfigsr*(j+1)+1)
                   
                    ## range *** when there are exactly 4 rows
                    range4=range(0, nfigsr*2+1) if j < 2 else range(nfigsr*2, n_figs+1)
            
                    ## set maxlistTT
                    maxlistTT=np.median(maxlistT)*0.3 if n_rows < 3 else min(maxlistT[min(range3):max(range3)]) if n_rows == 3 else min(maxlistT[min(range4):max(range4)])

                    ## If maxlistTT is less than then double
                    maxlistTT=maxlistTT*0.5 if maxlistTT < 300 and n_rows > 2 else maxlistTT

                    ## maxlistTT only for HTA1-HTA2 -- UGLY
                    maxlistTT=maxlistTT*0.55 if maxlistTT < 200 and n_rows > 2 and meta_of_step["pairs display"][0] == "HTA1-HTA2" else maxlistTT*0.2 if maxlistTT > 1000 and n_rows > 2 and meta_of_step["pairs display"][0] == "HTA1-HTA2" else maxlistTT
                
                    ## number of figs per "label" --> for the first unique label 
                    n_cols=len(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]])
                    
                    locals()["ax" + str(j)] = subplots[j].subplots(1, n_cols, sharex = True, sharey=True, squeeze = False)
                    locals()["ax" + str(j)] = locals()["ax" + str(j)].ravel()
                    
                    for i in range(0, n_cols):
                        test=imread(fname = meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['gfp path'].reset_index(drop = True)[i])
                        locals()["ax" + str(j)][i].imshow(X=test
                                      , cmap = mycmap
                                      , vmax = maxlistTT
                                      )
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['label display'].reset_index(drop = True)[i], size=8) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].set_xticks([])
                        locals()["ax" + str(j)][i].set_yticks([])
                        
                fig.suptitle('   ', size=50)
                return fig
            
            else:
                j=0
                
                ## number of figs per "label" --> for the first unique label 
                n_cols=len(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]])##CHANGE HERE##
                
                fig, axs = plt.subplots(1, n_cols, sharex = True, sharey=True, squeeze = False, figsize=(5, 20))
                axs = axs.ravel()
                
                ## Loop per axs
                for i in range(0, n_cols):
                    test=imread(fname = meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['gfp path'].reset_index(drop = True)[i])
                    axs[i].imshow(X=test, cmap = mycmap
                                  , vmax = np.median(maxlistT)*0.3
                                  )
                    axs[i].set_title(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['label display'].reset_index(drop = True)[i], size=8)
                    axs[i].set_xticks([])
                    axs[i].set_yticks([])
                        
                fig.suptitle('   ', size=50)
                return fig

        else:
          
            ## Notification "figure" for no image replicate
            fig, ax = plt.subplots()
            ax.axis([0, 1, 0, 1])
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.spines['bottom'].set_color('white')
            ax.spines['top'].set_color('white') 
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.set_title('No Replicate 2', size=15)## CHANGE HERE ##

            return fig
        
    ## Set plot for replicate 3 ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def rep3():## CHANGE HERE ##
        meta_of_step=meta2()
        
        if len(meta_of_step[meta_of_step["replicate"] == "replicate3"]) > 0:## CHANGE HERE ##
            
            ## set meta for replicate 1
            meta_of_step=meta_of_step[meta_of_step["replicate"] == "replicate3"].sort_values(by='label', ascending=False).reset_index(drop = True)## CHANGE HERE ## 
            
            ## set gfp path for range
            gfp_path=meta_of_step["gfp path"]

            ## Get range of both labels
            ima_range=range(meta_of_step.index.min()
                          , meta_of_step.index.max() + 1
                          )
                          
            ## Get maxs limits
            maxlist=maxs3()
            for i in ima_range:
                im1=imread(fname = gfp_path[i])
                maxi=im1.max()
                maxlist.append(maxi)
                maxlistT=maxlist[1:]
            
            ## T n figures
            n_figs=max(ima_range)+1
            
            ## Set number of rows based on # of labels 
            n_rows=len(meta_of_step["label"].unique())
            
            if n_rows > 2 :
            
                ## Set fig --> for entire set of Axs
                fig = plt.figure(figsize=(10,10), layout = 'tight')
                
                ## Set number of rows based on # of labels 
                subplots = fig.subfigures(n_rows, 1)
                
                #### Set Loop per number of labels
                for j in range(0, n_rows):
                  
                    ## n figures per row
                    nfigsr=round(n_figs/n_rows)
        
                    ## range of nfigsr  --> n  *** when there are exactly 3 rows
                    if len(meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]]) > 4 and j < 2 :
                        range3=range(0,nfigsr*2+1)
                    elif len(meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]]) < 5 and j == 0 :
                        range3=range(0, nfigsr+1)
                    else : 
                        range3=range(nfigsr,n_figs+1)        
                    # range3=range(nfigsr*j,nfigsr*(j+1)+1)
                    
                    ## range *** when there are exactly 4 rows
                    range4=range(0, nfigsr*2+1) if j < 2 else range(nfigsr*2, n_figs+1)
            
                    ## set maxlistTT
                    maxlistTT=np.median(maxlistT)*0.3 if n_rows < 3 else min(maxlistT[min(range3):max(range3)]) if n_rows == 3 else min(maxlistT[min(range4):max(range4)])

                    ## If maxlistTT is less than then double
                    maxlistTT=maxlistTT*0.5 if maxlistTT < 300 and n_rows > 2 else maxlistTT

                    ## maxlistTT only for HTA1-HTA2 -- UGLY
                    maxlistTT=maxlistTT*0.55 if maxlistTT < 200 and n_rows > 2 and meta_of_step["pairs display"][0] == "HTA1-HTA2" else maxlistTT*0.2 if maxlistTT > 1000 and n_rows > 2 and meta_of_step["pairs display"][0] == "HTA1-HTA2" else maxlistTT
                
                    ## number of figs per "label" --> for the first unique label 
                    n_cols=len(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]])
                    
                    locals()["ax" + str(j)] = subplots[j].subplots(1, n_cols, sharex = True, sharey=True, squeeze = False)
                    locals()["ax" + str(j)] = locals()["ax" + str(j)].ravel()
                    
                    for i in range(0, n_cols):
                        test=imread(fname = meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['gfp path'].reset_index(drop = True)[i])
                        locals()["ax" + str(j)][i].imshow(X=test
                                      , cmap = mycmap
                                      , vmax = maxlistTT
                                      )
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['label display'].reset_index(drop = True)[i], size=8) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].set_xticks([])
                        locals()["ax" + str(j)][i].set_yticks([])
                        
                fig.suptitle('   ', size=50)
                return fig
            
            else:
                j=0
                
                ## number of figs per "label" --> for the first unique label 
                n_cols=len(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]])##CHANGE HERE##
                
                fig, axs = plt.subplots(1, n_cols, sharex = True, sharey=True, squeeze = False, figsize=(5, 20))
                axs = axs.ravel()
                
                ## Loop per axs
                for i in range(0, n_cols):
                    test=imread(fname = meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['gfp path'].reset_index(drop = True)[i])
                    axs[i].imshow(X=test, cmap = mycmap
                                  , vmax = np.median(maxlistT)*0.3
                                  )
                    axs[i].set_title(meta_of_step[meta_of_step["label"] == meta_of_step["label"].unique()[j]]['label display'].reset_index(drop = True)[i], size=8)
                    axs[i].set_xticks([])
                    axs[i].set_yticks([])
                        
                fig.suptitle('   ', size=50)
                return fig

        else:
          
            ## Notification "figure" for no image replicate
            fig, ax = plt.subplots()
            ax.axis([0, 1, 0, 1])
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.spines['bottom'].set_color('white')
            ax.spines['top'].set_color('white') 
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.set_title('No Replicate 3', size=15)## CHANGE HERE ##

            return fig

    ## Acknowledgement text
    @render.text
    def acknowledge():## CHANGE HERE ##
        return f"All the images for each paralog pair, per replicate, have been thresholded separately for visualization only. \n\n" + \
                f"Supplementary data files are available from here: \n\n" + \
                f"Rohan Dandage, Mikhail Papkov, Brittany M. Greco, Dmytro Fishman, Helena Friesen, Kyle Wang, \n" + \
                f"Erin Styles, Oren Kraus, Benjamin Grys, Charles Boone, Brenda Andrews,Leopold Parts, Elena Kuzmin" + \
                f" \n'Single-cell imaging of protein dynamics of paralogs reveals mechanisms of gene retention.'" + \
                f" \nbioRxiv (2023): 2023-11.  doi: https://doi.org/10.1101/2023.11.23.568466"


# Close app ####################################################################
app = App(app_ui, server)
