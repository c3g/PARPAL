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
from htmltools import HTML, div

## root of data  
DATAROOT=os.getenv("PARPAL_DATA_LOCATION", "/data")

## Load data and compute static values
meta=pd.read_table(f"{DATAROOT}/metadata.tsv")## CHANGE HERE ##

## add DATAROOT env var to gfp path.   
meta["gfp path"]= f"{DATAROOT}/" + meta['gfp path']

## Select genes of interes to choose from 
genesofi=sorted(meta["pairs display"].unique().tolist())

## Select my picture colors
mycmap=colors.LinearSegmentedColormap.from_list("custom", colors=["#262626",'#83f52c'], N=50)


################################################################################
##################### ASSUMPTIONS FOR THIS TO WORK #############################

# 0 - Metadata must be in /data folder, along with all images and tables.  
# 1 - There must always be two pair in each paralos pair.  
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
        ui.input_action_button("submit", "Submit", class_="btn-success").add_style("background:#7393B3; border:white; font-family:Helvetica Neue;"), 
        # ui.output_ui("compute"), ## loading button slowes down processes
        ui.output_text_verbatim("C3G", placeholder = False).add_style("text-align:center; background:#F8F8F8; border:white; font-size:9pt; font-family:Helvetica Neue;"),
        ui.h3(
          HTML("<a href='https://computationalgenomics.ca/team_profiles/#GZ'>https://computationalgenomics.ca</a>")
            ).add_style("text-align:center; font-size:8pt; font-family:Helvetica Neue;"),  
        ui.output_image("C3Gimage", inline = True).add_style("text-align:center; padding:10px;"),
        shinyswatch.theme.litera(),
        width=400,
        padding=40, ## Padding btw things inside
    ),
    ## Beginning of the main page - Intro
    ui.output_image("EKimage", inline = True).add_style("text-align:center; padding-right:100px; padding-left:100px;"),
    ui.h3("Introduction").add_style("text-align:center; padding:40px; font-family:Helvetica Neue; font-size:22pt;"),
    ui.div(), 
    ui.output_text_verbatim("message", placeholder = True
        ).add_style("text-align:center; padding:10px; background:white; font-size:14pt; border:white; font-family:Helvetica Neue;"),
    ui.h6(
      "Email: ", HTML("<a href='mailto:elena.kuzmin@concordia.ca'>elena.kuzmin@concordia.ca</a>")
        ).add_style("text-align:center; font-size:11pt;"),
    ui.h6(  
      "Website: ", 
      HTML("<a href='https://kuzmin-lab.github.io/'>https://kuzmin-lab.github.io/</a>") 
        ).add_style("text-align:center; font-size:11pt; font-family:Helvetica Neue;"),
    ui.h1(" ").add_style("text-align:center; padding-top:100px; font-family:Helvetica Neue;"),
    ui.h6("Please allow for some loading time when switching between paralog pairs."
        ).add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue;"),
    ## Beginning of the tabs
    ui.page_navbar(
        ui.nav_panel("Scores", 
            ui.output_table("paralogredis"
            ).add_style("text-align:center; padding-top:30px; padding-bottom:50px; padding-right:300px; padding-left:300px; font-size:9pt; font-family:Helvetica Neue;"),
        ),
        ui.nav_panel("Paralog 1",
            ui.h6("Wild-type background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:0px;"),
            ui.output_plot("pair2_1", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:0px;"),
            ui.hr(), 
            ui.h6("Deletion background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:0px;"), 
            ui.output_plot("pair2_2", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:0px;"),
        ),
        ui.nav_panel("Paralog 2",
            ui.h6("Wild-type background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:0px;"),
            ui.output_plot("pair1_1", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:0px;"),
            ui.hr(), 
            ui.h6("Deletion background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:0px;"), 
            ui.output_plot("pair1_2", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:0px;"),
        ),
        ui.nav_panel("Supplemental files",
            ui.h5("Download below supp. data").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:50px;"),
            ui.div(),
            ui.download_link("T1", "Table S1 - Yeast strains and plasmids used in this study").add_style("font-size:11pt; font-family:Helvetica Neue;"),
            ui.div(),
            ui.download_link("T2", "Table S2 - Protein abundance").add_style("font-size:11pt; font-family:Helvetica Neue;"),
            ui.div(),
            ui.download_link("T3", "Table S3 - Manual inspection of paralog pairs and negative controls").add_style("font-size:11pt; font-family:Helvetica Neue;"),
            ui.div(),
            ui.download_link("T4", "Table S4 - The redistribution, relative abundance change and relocalization").add_style("font-size:11pt; font-family:Helvetica Neue;"),
            ui.div(),
            ui.download_link("T5", "Table S5 - Features of paralogs").add_style("font-size:11pt; font-family:Helvetica Neue;"),
            ui.div(),
            ui.download_link("D1", "Data  S1 - Features of single cells").add_style("font-size:11pt; font-family:Helvetica Neue;"),             
            # ui.h6("-- Features of single cells").add_style("font-size:10pt; font-family:Helvetica Neue;"),
            ui.div(""),
            ui.download_link("D2", "Data  S2 - Protein abundance per single cell").add_style("font-size:11pt; font-family:Helvetica Neue;"),
            ui.div(),
        ),
        id="tab",
    ).add_style("font-family:Helvetica Neue;"),
    ## Beginning of the closing statement. 
    ui.h1(" ").add_style("text-align:center; padding-top:100px; font-family:Helvetica Neue;"),
    ui.h3("References").add_style("text-align:center; padding:40px; font-family:Helvetica Neue;font-size:22pt;"),
    ui.output_text_verbatim("acknowledge"
        ).add_style("text-align:center; padding:10px; background:white; font-size:14pt; border:white; font-family:Helvetica Neue;"),
    ui.h6("DOI: ",
          HTML("<a href='https://doi.org/10.1101/2023.11.23.568466'>https://doi.org/10.1101/2023.11.23.568466</a>"),
        ).add_style("text-align:center; font-size:11pt; padding-bottom:50px; font-family:Helvetica Neue;"),
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
        return f"\n\n\n\n\n\n\n\n\n\n" + \
                f"Website developed by: \n" + \
                f"Gerardo Zapata, Rohan Dandage, \n" + \
                f"Vanessa Pereira and Elena Kuzmin \n\n" + \
                f"In Collaboraton with: \n" + \
                f"Canadian Centre for Computational Genomics (C3G) \n" 

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
                f"For details on the PARPAL project, data and website please contact Elena Kuzmin: \n"

    
    ## New mix table 
    @render.table
    @reactive.event(input.submit)
    def paralogredis():
        gene_of_step1=gene1()
        gene_of_step2=gene2()
        
        ## Set table options for 3 decimal points 
        pd.options.display.float_format = "{:,.3f}".format 
        
        redis=pd.read_csv(f"{DATAROOT}/scores/NewScores_forGerardo.tsv", sep='\t', float_precision='round_trip')
        redis=redis[(redis['Gene'] == gene_of_step1) | (redis['Gene'] == gene_of_step2)]
        redis["No redistribution or protein abundance change"] = ""
        
        redis=redis[['Paralog pair', 'ORF1-ORF2', 'Gene', 'ORF', 'LFC', 'q-value']] if len(redis) > 0 else redis[["No redistribution or protein abundance change"]]
        
        return redis 


    ## Set variable for max intensity value per replicate 
    @reactive.calc
    @reactive.event(input.submit)
    def maxs1_1():
        maxs=[0]
        return maxs

    ## Set variable for max intensity value per replicate 
    @reactive.calc
    @reactive.event(input.submit)
    def maxs1_2():
        maxs=[0]
        return maxs

    ## Set variable for max intensity value per replicate 
    @reactive.calc
    @reactive.event(input.submit)
    def maxs2_1():
        maxs=[0]
        return maxs

    ## Set variable for max intensity value per replicate 
    @reactive.calc
    @reactive.event(input.submit)
    def maxs2_2():
        maxs=[0]
        return maxs

    ## Set pair1 WT ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def pair1_1():## CHANGE HERE ##
        meta_of_step=meta2()
        meta_of_step=meta_of_step.reset_index(drop = True)
        
        ## Select First. ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]].sort_values(by='label', ascending=True).reset_index(drop = True) ## CHANGE HERE 
        
        ## Set inputs
        gfp_path=meta_of_step["gfp path"]
        
        ## Get range of both labels
        ima_range=range(meta_of_step.index.min()
                      , meta_of_step.index.max() + 1
                      )
        
        ## Get mins and maxs.  ## CHANGE HERE ##
        maxlist=maxs1_1()
        for i in ima_range:
            im1=imread(fname = gfp_path[i])
            maxi=np.percentile(im1, 99.9999)
            # maxi=im1.max()
            maxlist.append(maxi)
            maxlistT=maxlist[1:]
    
        ## make meta_of_step either WT or DELTA
        meta_of_step=meta_of_step[meta_of_step["status partner"] == "WT"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        # meta_of_step=meta_of_step[meta_of_step["status partner"] == "DELTA"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        
        ## get name of partner
        name=meta_of_step['label_display2'][0]
        
        ## Set number of rows/replicates within each partner 
        n_rows=len(meta_of_step['replicate'].unique())
        
        ###### ensure that there is at least one replicate
        if n_rows > 1:
          
            ## Set fig --> for entire set of Axs
            fig = plt.figure(constrained_layout=True,figsize=(10,10))
            
            ## Set number of rows based on # of reps
            subplots = fig.subfigures(n_rows, 1)
            
            #### Set Loop for each rep
            for j in range(0, n_rows):
              
                ## num  -- per rep
                number=int(float(j)+1)
                
                ## number of figs  -- per rep 
                n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
                
                #### if there is more than one image per replicate if statement
                if n_cols > 0:
                
                    locals()["ax" + str(j)] = subplots[j].subplots(1, n_cols, sharex = True, sharey=True, squeeze = False) if n_cols == 1 else subplots[j].subplots(1, n_cols, sharex = True, sharey=True)
                    locals()["ax" + str(j)] = locals()["ax" + str(j)].ravel() 
                
                    ## Set second loop for figs  -- per rep    
                    for i in range(0, n_cols):
                        test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].imshow(X=test, cmap = mycmap
                                      , vmax = np.percentile(maxlistT, 2.5)
                                      )
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].set_xticks([])
                        locals()["ax" + str(j)][i].set_yticks([])
                        locals()["ax" + str(j)][i].spines['bottom'].set_color('white')
                        locals()["ax" + str(j)][i].spines['top'].set_color('white')
                        locals()["ax" + str(j)][i].spines['right'].set_color('white')
                        locals()["ax" + str(j)][i].spines['left'].set_color('white')
                        
                
                #### else print message when there is no image within one replicate
                else: 
                    
                    locals()["ax" + str(j)] = plt.subplots()
                    locals()["ax" + str(j)].axis([0, 1, 0, 1])
                    locals()["ax" + str(j)].tick_params(axis='x', colors='white')
                    locals()["ax" + str(j)].tick_params(axis='y', colors='white')
                    locals()["ax" + str(j)].spines['bottom'].set_color('white')
                    locals()["ax" + str(j)].spines['top'].set_color('white') 
                    locals()["ax" + str(j)].spines['right'].set_color('white')
                    locals()["ax" + str(j)].spines['left'].set_color('white')
                    locals()["ax" + str(j)].set_title('No images for rep' + str(number) , size=15)## CHANGE HERE ##
        
                    
            return fig
            
        ## This is if there is one one rep  
        elif n_rows == 1:
            
            ## make new J    
            j=0
            
            ## num  -- per rep
            number=int(float(j)+1)
                
            ## number of figs  -- per rep 
            n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
            
            fig, axs = plt.subplots(1, n_cols, sharex = True, sharey=True, figsize=(5, 5))
            axs = axs.ravel()
            
            ## Loop per axs
            for i in range(0, n_cols):
                test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                axs[i].imshow(X=test, cmap = mycmap
                              , vmax = np.percentile(maxlistT, 2.5)
                              )
                axs[i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                axs[i].set_xticks([])
                axs[i].set_yticks([])
                axs[i].spines['bottom'].set_color('white')
                axs[i].spines['top'].set_color('white')
                axs[i].spines['right'].set_color('white')
                axs[i].spines['left'].set_color('white')
        
            return fig
          
        ###### else print message for all missing replicates - no figures
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
            ax.set_title('No replicates for ' + str(name) , size=15)## CHANGE HERE ##
            
            return fig

    ## Set pair1 DELTA ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def pair1_2():## CHANGE HERE ##
        meta_of_step=meta2()
        meta_of_step=meta_of_step.reset_index(drop = True)

        ## Select First. ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[0]].sort_values(by='label', ascending=True).reset_index(drop = True) ## CHANGE HERE 
        
        ## Set inputs
        gfp_path=meta_of_step["gfp path"]
        
        ## Get range of both labels
        ima_range=range(meta_of_step.index.min()
                      , meta_of_step.index.max() + 1
                      )
                      
        ## Get mins and maxs. ## CHANGE HERE ##
        maxlist=maxs1_2()
        for i in ima_range:
            im1=imread(fname = gfp_path[i])
            maxi=np.percentile(im1, 99.9999)
            maxlist.append(maxi)
            maxlistT=maxlist[1:]
        
        ## make meta_of_step either WT or DELTA. ## CHANGE HERE ##
        # meta_of_step=meta_of_step[meta_of_step["status partner"] == "WT"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        meta_of_step=meta_of_step[meta_of_step["status partner"] == "DELTA"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        
        ## get name of partner
        name=meta_of_step['label_display2'][0]
        
        ## Set number of rows/replicates within each partner 
        n_rows=len(meta_of_step['replicate'].unique())
        
        ###### ensure that there is at least one replicate
        if n_rows > 1:
            
            ## Set fig --> for entire set of Axs
            fig = plt.figure(constrained_layout=True,figsize=(10,10))
            
            ## Set number of rows based on # of reps
            subplots = fig.subfigures(n_rows, 1)
            
            #### Set Loop for each rep
            for j in range(0, n_rows):
              
                ## num  -- per rep
                number=int(float(j)+1)
                
                ## number of figs  -- per rep 
                n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
                
                #### if there is more than one image per replicate if statement
                if n_cols > 0:
                
                    locals()["ax" + str(j)] = subplots[j].subplots(1, n_cols, sharex = True, sharey=True, squeeze = False) if n_cols == 1 else subplots[j].subplots(1, n_cols, sharex = True, sharey=True)
                    locals()["ax" + str(j)] = locals()["ax" + str(j)].ravel() 
                
                    ## Set second loop for figs  -- per rep    
                    for i in range(0, n_cols):
                        test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].imshow(X=test, cmap = mycmap
                                      , vmax = np.percentile(maxlistT, 2.5)
                                      )
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].set_xticks([])
                        locals()["ax" + str(j)][i].set_yticks([])
                        locals()["ax" + str(j)][i].spines['bottom'].set_color('white')
                        locals()["ax" + str(j)][i].spines['top'].set_color('white')
                        locals()["ax" + str(j)][i].spines['right'].set_color('white')
                        locals()["ax" + str(j)][i].spines['left'].set_color('white')
                        
                
                #### else print message when there is no image within one replicate
                else: 
                    
                    locals()["ax" + str(j)] = plt.subplots()
                    locals()["ax" + str(j)].axis([0, 1, 0, 1])
                    locals()["ax" + str(j)].tick_params(axis='x', colors='white')
                    locals()["ax" + str(j)].tick_params(axis='y', colors='white')
                    locals()["ax" + str(j)].spines['bottom'].set_color('white')
                    locals()["ax" + str(j)].spines['top'].set_color('white') 
                    locals()["ax" + str(j)].spines['right'].set_color('white')
                    locals()["ax" + str(j)].spines['left'].set_color('white')
                    locals()["ax" + str(j)].set_title('No images for rep' + str(number) , size=15)## CHANGE HERE ##
        
                    
            return fig
            
        ## This is if there is one one rep  
        elif n_rows == 1:
            
            ## make new J    
            j=0
            
            ## num  -- per rep
            number=int(float(j)+1)
           
            ## number of figs  -- per rep 
            n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
            
            fig, axs = plt.subplots(1, n_cols, sharex = True, sharey=True, figsize=(5, 5))
            axs = axs.ravel()
            
            ## Loop per axs
            for i in range(0, n_cols):
                test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                axs[i].imshow(X=test, cmap = mycmap
                              , vmax = np.percentile(maxlistT, 2.5)
                              )
                axs[i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                axs[i].set_xticks([])
                axs[i].set_yticks([])
                axs[i].spines['bottom'].set_color('white')
                axs[i].spines['top'].set_color('white')
                axs[i].spines['right'].set_color('white')
                axs[i].spines['left'].set_color('white')
        
            return fig
          
        ###### else print message for all missing replicates - no figures
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
            ax.set_title('No replicates for ' + str(name) , size=15)## CHANGE HERE ##
            
            return fig


    ## Set pair2 WT ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def pair2_1():## CHANGE HERE ##
        meta_of_step=meta2()
        meta_of_step=meta_of_step.reset_index(drop = True)
        
        ## Select First ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[1]].sort_values(by='label', ascending=True).reset_index(drop = True) ## CHANGE HERE 
        
        ## Set inputs
        gfp_path=meta_of_step["gfp path"]
        
        ## Get range of both labels
        ima_range=range(meta_of_step.index.min()
                      , meta_of_step.index.max() + 1
                      )
                      
        ## Get mins and maxs. ## CHANGE HERE ##
        maxlist=maxs2_1()
        for i in ima_range:
            im1=imread(fname = gfp_path[i])
            maxi=np.percentile(im1, 99.9999)
            maxlist.append(maxi)
            maxlistT=maxlist[1:]
        
        ## make meta_of_step either WT or DELTA
        meta_of_step=meta_of_step[meta_of_step["status partner"] == "WT"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        # meta_of_step=meta_of_step[meta_of_step["status partner"] == "DELTA"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        
        ## get name of partner
        name=meta_of_step['label_display2'][0]
        
        ## Set number of rows/replicates within each partner 
        n_rows=len(meta_of_step['replicate'].unique())
        
        ###### ensure that there is at least one replicate
        if n_rows > 1:
            
            ## Set fig --> for entire set of Axs
            fig = plt.figure(constrained_layout=True,figsize=(10,10))
            
            ## Set number of rows based on # of reps
            subplots = fig.subfigures(n_rows, 1)
            
            #### Set Loop for each rep
            for j in range(0, n_rows):
              
                ## num  -- per rep
                number=int(float(j)+1)
                
                ## number of figs  -- per rep 
                n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
                
                #### if there is more than one image per replicate if statement
                if n_cols > 0:
                
                    locals()["ax" + str(j)] = subplots[j].subplots(1, n_cols, sharex = True, sharey=True, squeeze = False) if n_cols == 1 else subplots[j].subplots(1, n_cols, sharex = True, sharey=True)
                    locals()["ax" + str(j)] = locals()["ax" + str(j)].ravel() 
                
                    ## Set second loop for figs  -- per rep    
                    for i in range(0, n_cols):
                        test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].imshow(X=test, cmap = mycmap
                                      , vmax = np.percentile(maxlistT, 2.5)
                                      )
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].set_xticks([])
                        locals()["ax" + str(j)][i].set_yticks([])
                        locals()["ax" + str(j)][i].spines['bottom'].set_color('white')
                        locals()["ax" + str(j)][i].spines['top'].set_color('white')
                        locals()["ax" + str(j)][i].spines['right'].set_color('white')
                        locals()["ax" + str(j)][i].spines['left'].set_color('white')
                        
                
                #### else print message when there is no image within one replicate
                else: 
                    
                    locals()["ax" + str(j)] = plt.subplots()
                    locals()["ax" + str(j)].axis([0, 1, 0, 1])
                    locals()["ax" + str(j)].tick_params(axis='x', colors='white')
                    locals()["ax" + str(j)].tick_params(axis='y', colors='white')
                    locals()["ax" + str(j)].spines['bottom'].set_color('white')
                    locals()["ax" + str(j)].spines['top'].set_color('white') 
                    locals()["ax" + str(j)].spines['right'].set_color('white')
                    locals()["ax" + str(j)].spines['left'].set_color('white')
                    locals()["ax" + str(j)].set_title('No images for rep' + str(number) , size=15)## CHANGE HERE ##
        
                    
            return fig
            
        ## This is if there is one one rep  
        elif n_rows == 1:
            
            ## make new J    
            j=0

            ## num  -- per rep
            number=int(float(j)+1)
            
            ## number of figs  -- per rep 
            n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
            
            fig, axs = plt.subplots(1, n_cols, sharex = True, sharey=True, figsize=(5, 5))
            axs = axs.ravel()
            
            ## Loop per axs
            for i in range(0, n_cols):
                test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                axs[i].imshow(X=test, cmap = mycmap
                              , vmax = np.percentile(maxlistT, 2.5)
                              )
                axs[i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                axs[i].set_xticks([])
                axs[i].set_yticks([])
                axs[i].spines['bottom'].set_color('white')
                axs[i].spines['top'].set_color('white')
                axs[i].spines['right'].set_color('white')
                axs[i].spines['left'].set_color('white')
        
            return fig
          
        ###### else print message for all missing replicates - no figures
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
            ax.set_title('No replicates for ' + str(name) , size=15)## CHANGE HERE ##
            
            return fig

    ## Set pair2 DELTA ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def pair2_2():## CHANGE HERE ##
        meta_of_step=meta2()
        meta_of_step=meta_of_step.reset_index(drop = True)
        
        ## Select First. ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["Tpairs"] == meta_of_step["Tpairs"].unique()[1]].sort_values(by='label', ascending=True).reset_index(drop = True) ## CHANGE HERE 
        
        ## Set inputs
        gfp_path=meta_of_step["gfp path"]
        
        ## Get range of both labels
        ima_range=range(meta_of_step.index.min()
                      , meta_of_step.index.max() + 1
                      )
                      
        ## Get mins and maxs. ## CHANGE HERE ##
        maxlist=maxs2_2()
        for i in ima_range:
            im1=imread(fname = gfp_path[i])
            maxi=np.percentile(im1, 99.9999)
            maxlist.append(maxi)
            maxlistT=maxlist[1:]
        
        ## make meta_of_step either WT or DELTA. ## CHANGE HERE ##
        # meta_of_step=meta_of_step[meta_of_step["status partner"] == "WT"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        meta_of_step=meta_of_step[meta_of_step["status partner"] == "DELTA"].sort_values(by='image id', ascending=True).reset_index(drop = True)
        
        ## get name of partner
        name=meta_of_step['label_display2'][0]
        
        ## Set number of rows/replicates within each partner 
        n_rows=len(meta_of_step['replicate'].unique())
        
        ###### ensure that there is at least one replicate
        if n_rows > 1:
            
            ## Set fig --> for entire set of Axs
            fig = plt.figure(constrained_layout=True,figsize=(10,10))
            
            ## Set number of rows based on # of reps
            subplots = fig.subfigures(n_rows, 1)
            
            #### Set Loop for each rep
            for j in range(0, n_rows):
              
                ## num  -- per rep
                number=int(float(j)+1)
                
                ## number of figs  -- per rep 
                n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
                
                #### if there is more than one image per replicate if statement
                if n_cols > 0:
                
                    locals()["ax" + str(j)] = subplots[j].subplots(1, n_cols, sharex = True, sharey=True, squeeze = False) if n_cols == 1 else subplots[j].subplots(1, n_cols, sharex = True, sharey=True)
                    locals()["ax" + str(j)] = locals()["ax" + str(j)].ravel() 
                
                    ## Set second loop for figs  -- per rep    
                    for i in range(0, n_cols):
                        test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].imshow(X=test, cmap = mycmap
                                      , vmax = np.percentile(maxlistT, 2.5)
                                      )
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                        locals()["ax" + str(j)][i].set_xticks([])
                        locals()["ax" + str(j)][i].set_yticks([])
                        locals()["ax" + str(j)][i].spines['bottom'].set_color('white')
                        locals()["ax" + str(j)][i].spines['top'].set_color('white')
                        locals()["ax" + str(j)][i].spines['right'].set_color('white')
                        locals()["ax" + str(j)][i].spines['left'].set_color('white')
                        
                
                #### else print message when there is no image within one replicate
                else: 
                    
                    locals()["ax" + str(j)] = plt.subplots()
                    locals()["ax" + str(j)].axis([0, 1, 0, 1])
                    locals()["ax" + str(j)].tick_params(axis='x', colors='white')
                    locals()["ax" + str(j)].tick_params(axis='y', colors='white')
                    locals()["ax" + str(j)].spines['bottom'].set_color('white')
                    locals()["ax" + str(j)].spines['top'].set_color('white') 
                    locals()["ax" + str(j)].spines['right'].set_color('white')
                    locals()["ax" + str(j)].spines['left'].set_color('white')
                    locals()["ax" + str(j)].set_title('No images for rep' + str(number) , size=15)## CHANGE HERE ##
        
            return fig
            
        ## This is if there is one one rep  
        elif n_rows == 1:
            
            ## make new J    
            j=0

            ## num  -- per rep
            number=int(float(j)+1)
            
            ## number of figs  -- per rep 
            n_cols=len(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]])##CHANGE HERE##
            
            fig, axs = plt.subplots(1, n_cols, sharex = True, sharey=True, figsize=(5, 5))
            axs = axs.ravel()
            
            ## Loop per axs
            for i in range(0, n_cols):
                test=imread(fname = meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['gfp path'].reset_index(drop = True)[i]) ##CHANGE HERE##
                axs[i].imshow(X=test, cmap = mycmap
                              , vmax = np.percentile(maxlistT, 2.5)
                              )
                axs[i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8) ##CHANGE HERE##
                axs[i].set_xticks([])
                axs[i].set_yticks([])
                axs[i].spines['bottom'].set_color('white')
                axs[i].spines['top'].set_color('white')
                axs[i].spines['right'].set_color('white')
                axs[i].spines['left'].set_color('white')
        
            return fig
          
        ###### else print message for all missing replicates - no figures
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
            ax.set_title('No replicates for ' + str(name) , size=15)## CHANGE HERE ##
            
            return fig

    
    
    ## Download D1
    @render.download
    def D1():
        path=os.path.join(f"{DATAROOT}/to_download/DataS1.tsv")
        return path

    ## Download D2
    @render.download
    def D2():
        path=os.path.join(f"{DATAROOT}/to_download/DataS2.tsv")
        return path
    
    ## Download T1
    @render.download
    def T1():
        path=os.path.join(f"{DATAROOT}/to_download/TableS1.xlsx")
        return path

    ## Download T2
    @render.download
    def T2():
        path=os.path.join(f"{DATAROOT}/to_download/TableS2.xlsx")
        return path

    ## Download T3
    @render.download
    def T3():
        path=os.path.join(f"{DATAROOT}/to_download/TableS3.xlsx")
        return path

    ## Download T4
    @render.download
    def T4():
        path=os.path.join(f"{DATAROOT}/to_download/TableS4.xlsx")
        return path

    ## Download T5
    @render.download
    def T5():
        path=os.path.join(f"{DATAROOT}/to_download/TableS5.xlsx")
        return path



    ## Acknowledgement text
    @render.text
    def acknowledge():## CHANGE HERE ##
        return f"All the images for each paralog for both genetic backgrounds have been visualized using the same intensity settings. \n\n" + \
                f"Supplementary data files are available from here: \n\n" + \
                f"Rohan Dandage, Mikhail Papkov, Brittany M. Greco, Dmytro Fishman, Helena Friesen, Kyle Wang, \n" + \
                f"Erin Styles, Oren Kraus, Benjamin Grys, Charles Boone, Brenda Andrews, Leopold Parts, Elena Kuzmin" + \
                f" \n'Single-cell imaging of protein dynamics of paralogs reveals mechanisms of gene retention.'" + \
                f" \nbioRxiv (2023): 2023-11."


# Close app ####################################################################
app = App(app_ui, server)
