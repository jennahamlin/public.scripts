##single file the file must be called app.R in order to use R Shiny

# user interface for chromopaint: UI generates the webpage that a vistor sees
#it is really html, that is generated from R
library(dplyr)
library(shiny)

#setwd("~/Desktop/jenna/Shiny/chromopaint.app.v1")

##will need to concatenate the data. 
##cat AN1f.2.1.nearest.tsv DAV1e.1.nearest.tsv DAV1e.nearest.tsv  > combined.nearest.tsv
##remove duplicates lines like the first header row of each individual file
## awk '!seen[$0]++' combined.nearest.tsv > combined.nearest.final.tsv

chromdata<- read.csv("combined.nearest.final.csv",header=T)

#Define UI for application that draws  painted chromosomes
ui <- fluidPage(
  
  #App title
  titlePanel("Saccharomyces cerevisiae chromopainting") ,
  
  #Sidebar layout iwth input definitions
  sidebarLayout(
    
    #sidebar panel for inputs
    sidebarPanel(
      
      #input: select a strain
      selectInput('refer',h2('Strain'), choices =levels(chromdata$Strain)),
      
      #adding table tag to the sidebar
        tags$table(width ="100%", 
          tags$th("Donor clades:",  colspan="3", style="font-size: 25px"),
          tags$tr(  
            tags$td("Beer 1", style="font-size:20px; color:#ffcf40"),
            tags$td("Beer 2", style="font-size:20px;color:#FF6F69"),
            tags$td("Brazil 1", style="font-size:20px; color:#d896ff")),
          tags$tr(  
            tags$td("Cachaca 1", style="font-size:20px;color:#3366ff"),
            tags$td("Cachca 2", style="font-size:20px; color:#40e0d0"),
            tags$td("China", style="font-size:20px; color:gray")),
          tags$tr(
            tags$td("EU oak", style="font-size:20px;color:#a6dba0"),
            tags$td("Japan oak 1", style="font-size:20px;color:#5ddb3e"),
            tags$td("Japan oak 2", style="font-size:20px;color:#acdb3e")),
          tags$tr(
            tags$td("Malaysia", style="font-size:20px;color:#308CDF"),
            tags$td("North Africa", style="font-size:20px;color:#8B4513"),
            tags$td("NC oak", style="font-size:20px;color:#5aae61" )),
          tags$tr(
            tags$td("PA oak", style="font-size:20px;color:#005a32"),
            tags$td("Sake", style="font-size:20px;color:#e41a1c" ),
            tags$td("Tawain", style="font-size:20px;color:#87CEFA")),
          tags$tr(
            tags$td("West Africa", style="font-size:20px;color:orange"),
            tags$td("Wine", style="font-size:20px;color:#762a83"))),
      br(),br(),br(),#add space
      h4("Phylogenetic relationships for S. cerevisiae populations:"),
      tags$img(src='07.22.18.SCbackbone.beer1and2.brazil1.cachaca1and2.png',  width ="100%")
      
    ),
    

    #main panel for displayiing outputs
    mainPanel(
            h4(textOutput('selected_refer')),
              br(),#add space
              tableOutput("results"),
              br(),br(),br(),
              plotOutput("chromplot", height=800), #adjust the amount of space taken up by the plots with height
              br(),br(),br(),br(),
              h4('Brief methods'),
              p('The raw data used to generate the chromopainting plots was download from the', 
                a(href = 'https://www.ebi.ac.uk/ena', 'European Nucleotide Archive'),
                'and mapped to the sacCer3 reference genome (S288c) from',
                a(href = 'http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/', 'UCSC'),
                'using standard practices. I then itertively painted chromosomes with',
                a(href = 'https://github.com/bensassonlab/scripts/tree/master/fastatools', 'faChromopaint'), 
                'using a subset of strains and setting a maximum within clade divergence to 0.2%. 
                See the phylogeny for strains used as  poplation references. Note that white columns are regions
                of the genome in which divergence was greater than 0.2% and black columns are regions of the genome
                with too much missing data.'),
            br(),
            p('Questions? email me: jennahamlin (at) gmail (dot) com, or for a short summary about me see ', 
            a(href = 'https://jennahamlin.github.io','here'))
    )
  )
)


#server for chromopaint: set of instructions written in R for the server to follow
#tell the server what to do when the vistor wants to change and view
##first part of this function dynamically places the user selected strain in the header after Chromopainting of the selected strain
server<-function(input, output) {
  output$selected_refer <- renderText({ 
    paste("Chromopainting of the selected strain:")
  })

  output$results <- renderTable({
    filtered <-
      chromdata %>%
      filter(Strain ==input$refer)%>%
      slice(1:1)%>%
      select(Strain, Isolated.from, Geographic.origin)
    filtered
  })
  
  output$chromplot <-renderPlot({
    
    par(mfrow=c(8,2),mar=c(1,1,1,1),xpd=T, yaxt='n',bty="n")
    
    a <- as.data.frame(filter(chromdata, chromdata$Strain %in% input$refer))

    for (chr in c("chr01.mfa", "chr02.mfa", "chr03.mfa", "chr04.mfa", "chr05.mfa", "chr06.mfa", "chr07.mfa", "chr08.mfa", "chr09.mfa",
                  "chr10.mfa", "chr11.mfa", "chr12.mfa", "chr13.mfa", "chr14.mfa", "chr15.mfa", "chr16.mfa")){
      plot(c(0,max(a$winpos)),c(0,100),type="n",main=sub(".mfa","",chr),yaxt='n',xaxt='n', adj=0)
      rect(a$winpos[a$infile==chr]-20000,0,a$winpos[a$infile==chr],100,col=as.vector(a$cladecolor[a$infile==chr]))
    }
    
   })
  }

shinyApp(ui=ui, server = server)

##notes to myself. Have to publish via shinyapps.io account which requires the package 
#rsconnect. In order to publish and deploy the app you have to provide a token from 
#the shinyapps.io site. My shinyapps.io account is: https://jhamlin.shinyapps.io
#and the specific account for this app is https://jhamlin.shinyapps.io/chromopaintApp/