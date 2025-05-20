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
    # INTRO - MAIN PAGE #
    ui.output_image("EKimage", inline = True).add_style("text-align:center;padding-right:100px; padding-left:250px;"),
    # ui.h3("Introduction").add_style("text-align:center; padding:40px; font-family:Helvetica Neue; font-size:22pt;"),
    # ui.div(), 
    ui.output_text_verbatim("message", placeholder = True
        ).add_style("text-align:center; padding:5px; background:white; font-size:14pt; border:white; font-family:Helvetica Neue;"),
    ui.h6(" ").add_style("text-align:center; padding-top:20px; font-family:Helvetica Neue;"),
    ui.h6("Please allow for some loading time when switching between paralog pairs!"
        ).add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue;"),
    ui.h6("It is recommended to zoom in/out with cursor (mouse, trackpad) rather than browser (command +/-) "
        ).add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue;"),
    ui.h6("to avoid scale changes to the website."
        ).add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue;"),
    ## Beginning of the tabs
    ui.page_navbar(
        # GENE PAIR INFO - TAB 1 #
        ui.nav_panel("Scores", 
            ui.output_table("paralogredis", 
                ).add_style("padding-top:100px; padding-bottom:100px; font-family:Helvetica Neue; text-align:center;"),
            ui.div(),
            ui.output_text_verbatim("redisinfo", placeholder = True
                ).add_style("text-align:left; padding:10px; background:white; font-size:8pt; border:white; font-family:Helvetica Neue;"),
            ui.div(), 
        ),
        # GENE PAIR FIGURES WT - TAB 2 #
        ui.nav_panel("Paralog 1",
            ui.h6("Wild-type background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:0px;"),
            ui.output_plot("pair1_1", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:50px;"),
            ui.hr(), 
            ui.h6("Deletion background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:100px; padding-bottom:0px;"), 
            ui.output_plot("pair1_2", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:50px;"),
        ),
        # GENE PAIR FIGURES DELTA - TAB 3 #
        ui.nav_panel("Paralog 2",
            ui.h6("Wild-type background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:0px;"),
            ui.output_plot("pair2_1", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:50px;"),
            ui.hr(), 
            ui.h6("Deletion background").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:100px; padding-bottom:0px;"), 
            ui.output_plot("pair2_2", width = "1000px", height = "1000px").add_style("text-align:center; padding-top:50px;"),
        ),
        # EXTRA PAPER DOWNLOADS - TAB 4 #
        ui.nav_panel("Supplemental files",
            ui.h5("Click to download supp. data").add_style("text-align:center; font-size:12pt; font-family:Helvetica Neue; padding-top:50px; padding-bottom:50px;"),
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
        # OUTSIDE DATA - TAB 5 #
        ui.nav_panel("External data and resources", 
            ui.h6("SGD (Saccharomyces Genome Database): ").add_style("text-align:center; padding-top:100px; font-size:11pt; font-family:Helvetica Neue;"),
            ui.output_ui("YLink1").add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.div(), 
            ui.h6("Trigenic interaction fraction using single and double gene deletion mutants"
                ).add_style("padding-top:40px; text-align:center; font-size:11pt; font-family:Helvetica Neue;"),
            ui.output_ui("TI1").add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.output_ui("TIL").add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.div(),
            ui.h6("Protein-protein interaction change in response to paralog deletion"
                ).add_style("padding-top:40px; text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.output_ui("PP1").add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.output_ui("PPL").add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.div(),
            ui.h6("Regulation of protein level in response to paralog deletion"
                ).add_style("padding-top:40px; text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.output_ui("RP1").add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.output_ui("RPL").add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
            ui.div(),
        ),
        # REFERENCES - TAB 6 #
        ui.nav_panel("References",
            ui.h6(" ").add_style("text-align:center; padding-top:20px; font-family:Helvetica Neue;"),
            ui.output_text_verbatim("acknowledge1", placeholder = True
                ).add_style("text-align:center; background:white; font-size:12pt; padding-bottom:0px; border:white; font-family:Helvetica Neue;"),
            ui.h6("DOI: ", HTML("<a href='https://doi.org/10.1101/2023.11.23.568466'>https://doi.org/10.1101/2023.11.23.568466</a>"),
                ).add_style("text-align:center; font-size:8pt; padding-top:0px; font-family:Helvetica Neue;"),
            ui.h6(" ").add_style("text-align:center; padding-top:20px; font-family:Helvetica Neue;"),
            ui.output_text_verbatim("acknowledge2", placeholder = True
                ).add_style("text-align:center; background:white; font-size:12pt; padding-bottom:0px; border:white; font-family:Helvetica Neue;"),
            ui.h6("DOI: ", HTML("<a href='https://doi.org/10.1101/2025.03.04.641431'>https://doi.org/10.1101/2025.03.04.641431</a>"),
                ).add_style("text-align:center; font-size:8pt; padding-top:0px; padding-bottom:50px; font-family:Helvetica Neue;"),
            ui.h6(" ").add_style("text-align:center; padding-top:20px; font-family:Helvetica Neue;"),
            ui.output_text_verbatim("messagef", placeholder = True
                ).add_style("text-align:center; background:white; font-size:12pt; border:white; font-family:Helvetica Neue;"),
            ui.h6("Email: ", HTML("<a href='mailto:elena.kuzmin@concordia.ca'>elena.kuzmin@concordia.ca</a>")
                ).add_style("text-align:center; font-size:10pt;"),
            ui.h6("Website: ",HTML("<a href='https://kuzmin-lab.github.io/'>https://kuzmin-lab.github.io/</a>") 
                ).add_style("text-align:center; font-size:10pt; font-family:Helvetica Neue;"),
        ),
        id="tab",
    ).add_style("font-family:Helvetica Neue;"),
    # REFERENCES - MAIN PAGE #
    ## Beginning of the closing statement. 
    ui.h1(" ").add_style("text-align:center; padding-top:100px; font-family:Helvetica Neue;"),
    # ui.h3("References").add_style("text-align:center; padding:40px; font-family:Helvetica Neue;font-size:22pt;"),
    shinyswatch.theme.litera(),
    title="PARPAL",## Web title
)

# Server ######################################################################
def server(input, output, session):



    # ESSETNAIL VARIABLES#
    ## Add Image for side bar
    @render.image
    def C3Gimage():
        img: ImgData = {"src": f"{DATAROOT}/images/c3g.jpg", "width": "50%"}
        return img
    
    ## Set code block for sidebar
    @render.text
    def C3G():
        return f"\n\n\n\n\n" + \
                f"Website developed & maintained by: \n" + \
                f"Gerardo Zapata, Brittany Greco, \n" + \
                f"Rohan Dandage, Vanessa Pereira and Elena Kuzmin \n\n" + \
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
        img: ImgData = {"src": f"{DATAROOT}/images/20240417_PARPAL_logo.png", "width": "55%"}
        return img
    
    ## remake metadata table --> for gene pair
    @reactive.calc
    def meta2():
        meta2=meta[meta["pairs display"] == input.genepairs()].sort_values(by='label', ascending=False).reset_index(drop = True)
        return meta2

    ## make variable for the two genes 
    @reactive.calc
    @reactive.event(input.submit)
    def gene1():
        meta_gene=meta2()
        meta_gene[["label1", "label2"]]=meta_gene['label_display2'].str.split(" ", expand=True)
        lab1=meta_gene['label1'].str.contains("-GFP")
        
        gene1=meta_gene[lab1]["gene symbol query"].unique()[0]

        return gene1

    @reactive.calc
    @reactive.event(input.submit)
    def gene2():
        meta_gene=meta2()
        meta_gene[["label1", "label2"]]=meta_gene['label_display2'].str.split(" ", expand=True)
        lab2=meta_gene['label2'].str.contains("-GFP")
        
        gene2=meta_gene[lab2]["gene symbol query"].unique()[0]

        return gene2
      
    ## make variable for the two genes
    @reactive.calc
    @reactive.event(input.submit)
    def gene11():
        meta_gene=meta2()
        genes_of_step=meta_gene["pairs display"][1].split("-")
        gene11=genes_of_step[0]
        return gene11

    @reactive.calc
    @reactive.event(input.submit)
    def gene22():
        meta_gene=meta2()
        genes_of_step=meta_gene["pairs display"][1].split("-")
        gene22=genes_of_step[1]
        return gene22
    
    
    
    # INTRO - MAIN PAGE #
    ## Set code block for text example
    @render.text
    def message():
        return f"PARPAL is a web database for single-cell imaging \n" + \
                f"of protein dynamics of paralogs revealing sources of gene retention \n\n" + \
                f"Database statistics:  \n" + \
                f"Proteins screened = 164 \n" + \
                f"Paralog pairs screened = 82 \n" + \
                f"Total micrographs = ~3.5K \n" + \
                f"Total cells = ~460K \n"


    
    # GENE PAIR INFO - TAB 1 #
    ## New mix table 
    @render.table
    @reactive.event(input.submit)
    def paralogredis():
        
        redis=pd.read_csv(f"{DATAROOT}/scores/TableS4_forGerardo_qvalueedited.csv", sep=',')
        redis=redis[(redis['Gene'] == gene11()) | (redis['Gene'] == gene22())]
        redis["No redistribution or protein abundance change"] = ""
        redis=redis.rename({'Gene':'Gene-GFP'
                ,'ORF':'ORF-GFP'
                ,'redistribution score':'redistribution score*'
                ,'redistribution':'redistribution*'
                ,'relative abundance change, LFC':'relative abundance change, LFC**'
                ,'relative abundance change, q-value':'relative abundance change, q-value**'
                ,'relative abundance change type':'relative abundance change type**'
                , 'relocalization type':'relocalization type***'
                , 'relocalization description':'relocalization description***'}
            , axis='columns')
        # redis.loc[redis['Paralog pair']=='CUE1-CUE4', ['relative abundance change, q-value']] = 0
        
        ## Loop to fix notation for CUE1 or POR1
        if (gene11() == "CUE1" or gene11() == "POR1"):
            
            ## Set q-value column as float with 1 digit
            redis['relative abundance change, q-value**']=redis['relative abundance change, q-value**'].map('{:.0f}'.format)
            
        else:

            ## Set q-value column as scientific notation (3 sig digits) with 'x10' instead of 'e'
            redis['relative abundance change, q-value**']=redis['relative abundance change, q-value**'].map('{:.3e}'.format).replace("e", "x10^",regex=True)
        
        ## Loop to display whole table or only the one column
        if len(redis) > 0:
            
            ## Selec column
            redis=redis.drop("No redistribution or protein abundance change", axis=1)
          
            ## Change name - center
            redis=redis.style.format({'relative abundance change, LFC**': '{:.3f}',
                                      'redistribution score*': '{:.3f}'})\
                .set_properties(**{'text-align':'center', 'font-size':'10px', 
                                                'font-family':'Helvetica Neue',}) \
                .hide(axis='index')\
                .set_table_styles([{'selector':'th', 'props': 
                                          'text-align:center; font-size:12px; font-family:Helvetica Neue; border-bottom:1px solid grey; border-left:0.5px solid #eee; border-right:0.5px solid #eee; border-top:0.5px solid #eee'},
                                    {'selector':'th:hover', 'props': 'background-color: #ebfce8;'},
                                    {'selector':'td:hover', 'props': 'background-color: #ebfce8;'},
                                    {'selector': 'td', 'props': 'border-right: 0.5px solid #eee; border-left: 0.5px solid #eee; border-bottom: 0.5px solid #eee;'}])
                                    
            # ## Change name - center
            # redis=redis.style.format({'relative abundance change, LFC**': '{:.3f}',
            #                           'redistribution score*': '{:.3f}'})\
            #     .set_properties(**{'text-align':'center', 'font-size':'10px', 
            #                                     'font-family':'Helvetica Neue',}) \
            #     .hide(axis='index')\
            #     .set_table_styles([{'selector':'th', 'props': 
            #                               'text-align:center; font-size:12px; font-family:Helvetica Neue; border-bottom:2px solid black'},
            #                         {'selector': 'td', 'props': 'border-right: 1px solid black;'}])

        else:
            
            ## Selec column
            redis=redis[["No redistribution or protein abundance change"]]
            
            ## Change name - center
            redis=redis.style.set_properties(**{'text-align': 'center', 'font-size': '20px', 'font-family': 'Helvetica Neue'})\
                .set_table_styles([{'selector':'th', 'props': 'text-align: center; font-size: 20px; font-family:Helvetica Neue;'}])
        
        return redis
    
    # GENE PAIR INFO - TAB 1 #
    ## Supp info for redistribution table
    @render.text
    @reactive.event(input.submit)
    def redisinfo():
      
        redis=pd.read_csv(f"{DATAROOT}/scores/TableS4_forGerardo_qvalueedited.csv", sep=',')
        redis=redis[(redis['Gene'] == gene11()) | (redis['Gene'] == gene22())]
        redis["No redistribution or protein abundance change"] = ""
        
        ## If there is anything in redis then print, if not then print 'space'
        if len(redis) > 0:
        
            return f"* - Derived from the analysis of deep neural network features \n" + \
                    f"** - Derived from mean GFP pixel intensity \n" + \
                    f"*** - Derived from scoring by visual inspection \n" + \
                    f"See references tab for more details in the corresponding manuscripts."
        
        else:
          
            return f" "


    # GENE PAIR FIGURES - ESSENTIALS - tab 2#
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

    # GENE PAIR FIGURES WT - TAB 2 #
    ## Set pair1 WT ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def pair1_1():## CHANGE HERE ##
        meta_of_step=meta2()
        meta_of_step=meta_of_step.reset_index(drop = True)
        
        ## Select First gene pair. ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["gene symbol query"] == gene1()].sort_values(by=['label display'], ascending=False).reset_index(drop = True) ## CHANGE HERE
        
        ## LOOP only for NSG1-NSG2 - gene pari 1 - remove rep3
        if meta_of_step['pairs display'].unique() == "NSG1-NSG2":
            meta_of_step=meta_of_step[meta_of_step["replicate"] != 'replicate3'].sort_values(by=['label display'], ascending=False).reset_index(drop = True) ## CHANGE HERE
        else:
            meta_of_step=meta_of_step
       
        ## Set inputs
        if meta_of_step['pairs display'].unique() == "RPS22A-RPS22B":
            meta_gfp=meta2()
            gfp_path=meta_gfp["gfp path"]
        else:
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
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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
                axs[i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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

        ## Select First gene pair. ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["gene symbol query"] == gene1()].sort_values(by=['label display'], ascending=False).reset_index(drop = True) ## CHANGE HERE
      
        ## LOOP only for NSG1-NSG2 - gene pari 1 - remove rep3
        if meta_of_step['pairs display'].unique() == "NSG1-NSG2":
            meta_of_step=meta_of_step[meta_of_step["replicate"] != 'replicate3'].sort_values(by=['label display'], ascending=False).reset_index(drop = True) ## CHANGE HERE
        else:
            meta_of_step=meta_of_step
       
        ## Set inputs
        if meta_of_step['pairs display'].unique() == "RPS22A-RPS22B":
            meta_gfp=meta2()
            gfp_path=meta_gfp["gfp path"]
        else:
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
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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
                axs[i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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



    # GENE PAIR FIGURES DELTA - TAB 3 #
    ## Set pair2 WT ## CHANGE HERE ##
    @render.plot
    @reactive.event(input.submit)
    def pair2_1():## CHANGE HERE ##
        meta_of_step=meta2()
        meta_of_step=meta_of_step.reset_index(drop = True)
        
        ## Select SECOND gene pair. ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["gene symbol query"] == gene2()].sort_values(by=['label display'], ascending=False).reset_index(drop = True) ## CHANGE HERE
        
        ## Set inputs
        if meta_of_step['pairs display'].unique() == "RPS22A-RPS22B":
            meta_gfp=meta2()
            gfp_path=meta_gfp["gfp path"]
        else:
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
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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
                axs[i].set_title("\n\n\n" + meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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
        
        ## Select SECOND gene pair. ## CHANGE HERE ##
        meta_of_step=meta_of_step[meta_of_step["gene symbol query"] == gene2()].sort_values(by=['label display'], ascending=False).reset_index(drop = True) ## CHANGE HERE
          
        ## Set inputs
        if meta_of_step['pairs display'].unique() == "RPS22A-RPS22B":
            meta_gfp=meta2()
            gfp_path=meta_gfp["gfp path"]
        else:
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
                        locals()["ax" + str(j)][i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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
                axs[i].set_title(meta_of_step[meta_of_step["replicate"] == meta_of_step["replicate"].unique()[j]]['label display'].reset_index(drop = True)[i] + " - rep" + str(number), size=8, y=0, pad=-2, verticalalignment="top") ##CHANGE HERE##
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

    
    
    # EXTRA PAPER DOWNLOADS - TAB 4 #
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



    # OUTSIDE DATA - TAB 5 #
    ## Yeast genome Link 1
    # @reactive.calc
    # @reactive.event(input.submit)
    # def YL1():
    #     urlYL1="https://www.yeastgenome.org/locus/" + gene1()
    #     return ui.tags.a("Click here" , href=urlYL1, target='_blank')
    # 
    # ## actual sentence YLink1
    # @render.ui
    # @reactive.event(input.submit)
    # def YLink1():
    #     return ui.div(YL1(), " - ", gene11())

    ## Links to Sacc websites
    @render.ui
    @reactive.event(input.submit)
    def YLink1():
        
        ## URL to each sacc. site      
        urlYL1="https://www.yeastgenome.org/locus/" + gene1()
        urlYL2="https://www.yeastgenome.org/locus/" + gene2()
        
        YLone=ui.tags.a(gene11(), href=urlYL1, target="_blank")
        YLtwo=ui.tags.a(gene22(), href=urlYL2, target="_blank")
 
        return ui.div(YLone, HTML("<br>"), YLtwo)


    ### Transgenic interaction
    ## ORF1 or ORF2 Value
    @render.ui
    @reactive.event(input.submit)
    def TI1():

        ## Read
        TItable=pd.read_csv(f"{DATAROOT}/external/Kuzmin_et_al_2020.csv", sep=',')
        TItable=TItable[(TItable['Gene1'] == gene11()) | (TItable['Gene2'] == gene22()) | (TItable['Gene2'] == gene11()) | (TItable['Gene1'] == gene22())][['Trigenic interaction fraction class']].values.tolist()[0]

        return ui.div(ui.div("Trigenic interaction fraction class: ") , TItable)

    ## Link
    @render.ui
    def TIL():
        urlTIL1="https://www.science.org/doi/10.1126/science.aaz5667"
        urlTIL2="http://boonelab.ccbr.utoronto.ca/paralogs/"
        
        one=ui.tags.a("Kuzmin et al. Science 2020", href=urlTIL1, target="_blank")
        two=ui.tags.a("Trigenic interaction data portal", href=urlTIL2, target="_blank")
        
        return ui.div(one, " or ", two)
      
    ### Protein-Protein interaction
    ## ORF1 or ORF2 Conclusion
    @render.ui
    @reactive.event(input.submit)
    def PP1():

        ## Read
        PPtable=pd.read_csv(f"{DATAROOT}/external/Diss_et_al_2017_250123.csv", sep=',')
        PPtable=PPtable[(PPtable['Gene1'] == gene11()) | (PPtable['Gene2'] == gene22()) | (PPtable['Gene2'] == gene11()) | (PPtable['Gene1'] == gene22())][['Conclusion']].values.tolist()[0]

        return ui.div(ui.div("Conclusion: ") , PPtable)

    ## Link
    @render.ui
    def PPL():
        urlPPL="https://www.science.org/doi/full/10.1126/science.aai7685"
        
        return ui.tags.a("Diss et al Science 2017", href=urlPPL, target="_blank")
    
    ## Regulation Protein Level
    # ORF1 or ORF2 Value
    @render.ui
    @reactive.event(input.submit)
    def RP1():

        ## Read
        RPtable=pd.read_csv(f"{DATAROOT}/external/DeLuna_et_al_2010.csv", sep=',')
        RPtable=RPtable[(RPtable['Gene1'] == gene11()) | (RPtable['Gene2'] == gene22()) | (RPtable['Gene2'] == gene11()) | (RPtable['Gene1'] == gene22())][['Response']].values.tolist()[0]

        return ui.div(ui.div("Response: ") , RPtable)

    ## Link
    @render.ui
    def RPL():
        urlRPL="https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000347"
        
        return ui.tags.a("DeLuna et al PLoS Biology 2010", href=urlRPL, target="_blank")



    # REFERENCES - TAB 6 #
    ## Acknowledgement text1
    @render.text
    def acknowledge1():## CHANGE HERE ##
        return f"\nAll images for each paralog, for both genetic backgrounds have been visualized using the same intensity settings. \n\n\n\n" + \
                f"Supplementary data files are available from: \n\n" + \
                f"Rohan Dandage, Mikhail Papkov, Brittany M. Greco, Dmytro Fishman, Helena Friesen, Kyle Wang, \n" + \
                f"Erin Styles, Oren Kraus, Benjamin Grys, Charles Boone, Brenda Andrews, Leopold Parts, Elena Kuzmin" + \
                f" \n'Single-cell imaging of protein dynamics of paralogs reveals mechanisms of gene retention.'" + \
                f" \nbioRxiv (2023): 2023-11."
                
    ## Acknowledgement text2
    @render.text
    def acknowledge2():## CHANGE HERE ##
        return f"Brittany M. Greco*, Gerardo Zapata*, Rohan Dandage, Mikhail Papkov, \n" + \
                f"Vanessa Pereira, François Lefebvre, Guillaume Bourque, Leopold Parts and Elena Kuzmin" + \
                f" \n'PARPAL: PARalog Protein Redistribution using Abundance and Localization in Yeast Database'" + \
                f" \nbioRxiv (2025)."

    ## Kuzmin Info
    @render.text
    def messagef():
        return f"For details on the PARPAL project, data and website please contact Elena Kuzmin: \n"




# Close app ####################################################################
app = App(app_ui, server)
