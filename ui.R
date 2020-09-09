########piRNA User Interface####
###Amhed Vargas
###amhed.velazquez@kaust.edu.sa
###UI

#Load libraries
library(shiny)
library(shinythemes)
library(DT)

# Define User interface
shinyUI(
    fluidPage(
        tags$head(
            tags$link(rel="stylesheet",type = "text/css", href="bootstrap.min.css")
        ),
        ##Costum extra styles: single sliders background and title of navbar 
        tags$style(type = 'text/css', 
                   ".js-irs-none .irs-single, .js-irs-none .irs-bar-edge, .js-irs-none .irs-bar {
                          background: transparent;
                          border-top-color: transparent;
                          border-bottom-color: transparent;
                          border-left-color: transparent;
                          border-right-color: transparent}
               .navbar-default .navbar-brand:hover {color: #ffffff;}
               "),
        tags$style(type='text/css', '#AdvancedFragment {white-space: pre-wrap;}'),
        tags$style(type='text/css', '#SimpleFragment {white-space: pre-wrap;}'),
        tags$style(type='text/css', '#piBoxes .form-group {margin-bottom: 0px; margin-top: 0px;}'),
        #Main tab pages
        navbarPage(
            title=actionLink("link_to_tabpanel_title", HTML("<b>piRNAi</b>")),
            windowTitle="piRNAi app",
            id = "panels",
            
            tabPanel("Designer",
                     mainPanel(
                         h1("Design piRNAi fragments"),
                         HTML("<h4>for <i>C. elegans</i> and <i>C. briggsae</i> <a href=\"https://wormbase.org/about/wormbase_release_WS270\">WS270</a>.</h4>"),
                         tabsetPanel(
                            tabPanel("Simple",
                                     br(),
                         textAreaInput("geneinput", label = "Target gene", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                         actionButton("actiongenesearch", label = "Search gene"),
                         hr(),
                         uiOutput("DesignControls"),
                         hr(),
                         verbatimTextOutput("ErrorMessage"),
                         #tableOutput(otherPis),
                         htmlOutput("SelPiTabSummary"),
                         tableOutput('SelPiTab'),
                         verbatimTextOutput("SimpleFragment"),
                         uiOutput("downloadseq")
                         ),
                         tabPanel("Advanced",
                         h3("Make your own construct:"),
                    fluidRow(
                        ##Chose cluster
                        
                        radioButtons("clustercon", label = HTML("Select piRNA cluster
                                                     [<a href=\"\" onclick=\"$('#explain_cluster_adv').toggle(); return false;\">info</a>]
                                                     "),
                                     choices = list("21ur-1224 (6 sites)" = 1, "21ur-1692 (6 sites)" = 2, "21ur-8831 (7 sites)" = 3, "21ur-1610 (8 sites)" = 4), selected = 1, width='100%', inline= TRUE),
                        HTML("
                     <p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_cluster_adv\">
            We recommend to use the cluster 21ur-1224.
                                     </div></p>
                     "),
                        conditionalPanel(condition = "input.clustercon==1",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_1", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==2",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_2", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==3",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6"),
                        textAreaInput("piRNAseq7_3", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 7")),
                        br()
                        )),
                        conditionalPanel(condition = "input.clustercon==4",
                        column(4,
                               tags$div(id="piBoxes",
                                        class="my_class",
                        textAreaInput("piRNAseq1_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 1"),
                        textAreaInput("piRNAseq2_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 2"),
                        textAreaInput("piRNAseq3_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 3"),
                        textAreaInput("piRNAseq4_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 4"),
                        textAreaInput("piRNAseq5_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 5"),
                        textAreaInput("piRNAseq6_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 6"),
                        textAreaInput("piRNAseq7_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 7"),
                        textAreaInput("piRNAseq8_4", label="", rows=1, cols=21, resize="none", placeholder = "piRNA 8")),
                        br()
                        )),
                        column(8,
                        verbatimTextOutput("AdvancedFragment"),
                        uiOutput("downloadconstruct"))
                        ),
                        
                     fluidRow(
                         column(2,actionButton("actionclean", label = "Reset")),
                         column(2,actionButton("actionconstruct", label = "Generate piRNAi cluster"))),
                         verbatimTextOutput("AdvancedErrorMessage"),
                         hr(),
                         textAreaInput("Advancedgeneinput", label = "Target gene", value = "", resize="none", placeholder= "WormbaseID, transcript or gene name", rows=1),
                         actionButton("actionAdvsearch", label = "Search for specific piRNAs"),
                         hr(),
                         uiOutput("AdvDesignControls"),
                         hr(),
                         DT::dataTableOutput('AllPiTab')
                     )
            ))),
            ###Background
            tabPanel("Background",
                     mainPanel(
                         h3("A database and cloning system for a genome-wide piRNA interference collection"),
                         HTML("<p align=\"justify\">
                         In <i>C. elegans</i>, genome-wide tools based on RNA interference have been used to do “systems biology” screens. The collection was based on the observation that bacteria expressing a double-stranded RNA (dsRNA) can elicit an RNAi response in worms if they eat the bacteria (Timmons and Fire, 1998). This lead to the creation of a genome-wide collection of bacteria expressing dsRNA against most of <i>C. elegans</i> genes (Ahringer lab library). This has been used a lot in the field but has some limitations. In the germline, the phenotype from ingested dsRNA is considerably weaker than injected dsRNA.<br> 
<br>
Recently, our lab has developed methods to silence genes via piRNAs instead (Priyadarshini <i>et al., in preparation</i>). piRNAs are a class of small RNAs that are active in the germline. We only recently learned how the piRNAs actually recognize genes (Heng-Chi lab paper and Mello lab paper). 
Knowing the rules, mean that we can design piRNAs to target specific genes and so, that is what we have done here.

This app helps to "),
                         actionLink("link_to_tabpanel_design", "design"),
                         HTML(" piRNAi fragments easily. Alternatively, you can "),
                         actionLink("link_to_tabpanel_downloads", "download"),
                         HTML(" our designs and see them in a genome browser."),
                     HTML(" </p>")
                     )
            ),
            ###About
            tabPanel("Downloads",
                     mainPanel(
                         h3("Tracks"),
                         HTML("<p align=\"justify\">
                         <a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/tracks/piRNAi/Celegans.tar.gz\"><i>C. elegans</i> piRNA sequences</a><br>
                         <a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/tracks/piRNAi/Cbriggsae.tar.gz\"><i>C. briggsae</i> piRNA sequences</a><br>
                      </p>")
                     )
            ),
            ###About
            tabPanel("About",
                     mainPanel(
                         h3("The app"),
                         HTML("<p align=\"justify\">This website is generated via custom modified css/html code running in R via the shiny library.
                 <br>All the templates, libraries, and programs used to produce this site are under the MIT and GNU licenses.</p>"),
                         h3("The piRNAi algorithm"),
                         HTML("<p align=\"justify\">
                      Deviced by Christian Frøkjær-Jensen, tested experimentally by Monika Priyadarshini and implemented computationally by <a href=\"https://www.researchgate.net/profile/Amhed_Vargas_Velazquez\">Amhed Missael Vargas Velazquez</a> </p>"),
                         h3("The Laboratory of Synthetic Genome Biology"),
                         HTML("<p align=\"justify\">
                 The Laboratory of Synthetic Genome Biology is located in building 2 - level 3 (Ibn Al-Haytham – Above Spine) at King Abdullah University of Science and Technology (KAUST).
                 <br><i>Contact info</i>:<br>Christian-Froekjaer Jensen, Ph.D. 
                 <br>Assistant Professor of Bioscience
                 <br>Laboratory of Synthetic Genome Biology
                 <br>Email: <a href=\"mailto:cfjensen@kaust.edu.sa\">cfjensen@kaust.edu.sa</a>
                      </p>")
                     )
            )
        ),
        hr(),
        HTML("<a href=\"https://syngenbio.kaust.edu.sa/\">Syntetic genome biology laboratory @KAUST</a><br>"),
        HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
        HTML("<a href=\"mailto:cfjensen@kaust.edu.sa\">Contact us!</a>")
        
    )
)
