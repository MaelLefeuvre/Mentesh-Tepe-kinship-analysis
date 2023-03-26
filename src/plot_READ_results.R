#!/usr/bin/env Rscript

# -------------------------------------------------------------------------------------------------------------------- #
#                                                plot_READ_results.R                                                   #
# -------------------------------------------------------------------------------------------------------------------- #
# Author     : Mael Lefeuvre | Description: generate pretty and synthetic visualization for READ. pairs of interest    #
# Creation   : 2021-06-24    |              are highlighted while others are aggregated within a single boxplot.       #
# -------------------------------------------------------------------------------------------------------------------- #
# Dependencies: [CRAN] Xmisc, plotly, htmlwidgets                                                                      #
# -------------------------------------------------------------------------------------------------------------------- #
# Usage: Rscript ./plot_READ_results.R --results ./READ_results --meansP0 ./meansP0_AncientDNA_normalized              #
# -------------------------------------------------------------------------------------------------------------------- #


# ---------------------------------------------------- Functions ----------------------------------------------------- #
load_dependencies = function(packages){
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)){
    install.packages(packages[!installed_packages])
  }
  lapply(packages, library, character.only = TRUE)
}

# ---- Argument parser
parse_arguments= function(){
    parser = Xmisc::ArgumentParser$new()
    
    parser$add_usage("plot_READ_results.R [options]")
    parser$add_description('Generate synthetic visualization for the READ Kinship estimation method')
    
    parser$add_argument('--results', type='character', default="./READ_results", help='[Required] (Path): main output dataset of READ kinship estimator. Most probably named "READ_results"'
    )
    parser$add_argument('--meansP0', type='character', default="./meansP0_AncientDNA_normalized",
			help='[Required] (Path): secondary output dataset of READ kinship estimator. Most probably named "meansP0_AncientDNA_normalized"')
    parser$add_argument('--regex', type='character',
			help='[Required] (String): regular expression used to find the main individual pairs of interest. [Example]: "MT.*MT.*" will match any pair of individuals both belonging to the MT population. (e.g. : "MT23MT26", "MT8MT32")')
    parser$add_argument('--mainName', type='character', default="Main", help='[Optional] (String): Full name of the population of interest. Merely for plotting purposes.')
    parser$add_argument('--proxyName', type='character', default="Proxy", help='[Optional] (String): Full name of the surrogate, or "proxy" population. Merely for plotting purposes.')

    parser$add_argument('--format', type='character', default='png', help='[Optional] (String) Output file format ("png", "svg").')

    parser$add_argument('--h', type='logical', action='store_true', help='Print this help message.')
    parser$add_argument('--help', type='logical', action='store_true', help='Print this help message.')
 
    parser$helpme()
    optargs=parser$get_args()
    return(optargs)
}

# ------------------------------------------------------ Main -------------------------------------------------------- #
main = function(){
  # options
  options(viewer=NULL)

  # Load dependencies 
  dependencies = c("plotly", "pracma", "Xmisc", "htmlwidgets")
  load_dependencies(dependencies)

  # Parse command-line arguments
  out <- tryCatch(
    {
      optargs=parse_arguments()
    },
    error=function(cond){
      message(paste(cond,'\n'))
      return(1)
    }
  )  

  # Import data
  results = read.table(optargs$results, header=T, sep="\t")
  meansP0 = read.table(optargs$meansP0, header=T, sep=" ")

  #Extract main individuals
  main_pairs    = grepl(optargs$regex, meansP0$PairIndividuals)
  out = tryCatch(
    expr = {
      !any(main_pairs)
    },
    error = function(cond){
      message(paste0("Error: Provided regex pattern: '",optargs$regex,"' failed to catch any pair of individuals."))
    }
  )
  proxy_pairs = !main_pairs

  # plot
  main_only         = meansP0[main_pairs,]
  main_only         = main_only[order(main_only$Normalized2AlleleDifference),]
  main_positions    = seq(sum(main_pairs))/10
  proxy_positions = max(main_positions) + 2*min(main_positions)
  
  fig <- plot_ly()
  
  # Plot main individuals as separate markers, with error bars
  fig <- fig %>% add_markers(type    = 'scatter',
                             mode    = 'markers',
                             x       = main_positions, 
                             y       = main_only$Normalized2AlleleDifference,
                             name    = paste('<b>', optargs$mainName, 'individuals</b>'),
                             marker  = list(color = '#648FFF',
                                            size  = 12,
                                            line  = list(color='#785EF0')
                                           ),
                             error_y = ~list(array = main_only$StandardError*2,
                                             color = '#000000'
                                            ),
                           )
  
  # Add surrogate population individuals as a single boxplot
  fig <- fig %>% add_boxplot(x           = 1,
                             y           = meansP0$Normalized2AlleleDifference[proxy_pairs],
                             jitter      = 1,
                             pointpos    = 0,
                             boxpoints   = 'all',
                             marker      = list(color = 'rgba(100,100,100,0.3)'),
                             line        = list(color = 'rgba(7,40,89,0.6)'),
                             name        = paste0("<b>",optargs$proxyName," individuals</b>")
                            )
  
  # Create colored rectangles for each relatedness order.
  rectangles  = list()
  phi0        = c(0.90625, 0.8125, 0.625, 0)
  kinships    = c("Second Degree", "First Degree", "Identical Twins")
  colors      = c("#009E73","#F0E442","#CC79A7")
  
  for (i in 1:length(kinships)){
      rectangle=list()
      rectangle[["type"]]      = "rect"
      rectangle[["fillcolor"]] = colors[i]
      rectangle[["line"]]      = list(color=colors[i])
      rectangle[["opacity"]]   = 0.15
      rectangle[["x0"]]        = 0
      rectangle[["x1"]]        = 1
      rectangle[["xref"]]      = "paper"
      rectangle[["y0"]]        = phi0[i+1]
      rectangle[["y1"]]        = phi0[i]
      rectangles               = c(rectangles, list(rectangle))
  }
  
  
  # Add labels for each relatedness order
  fig <- fig %>% add_annotations(x         = 1,
                                 y         = phi0[1:length(phi0)-1],
                                 yanchor   = "top",
                                 xanchor   = "right",
                                 xref      = "paper",
                                 text      = sprintf("<b>%s</b>", kinships),
                                 showarrow = FALSE,
                                 font      = list(size=20)
                                )
  
  fig <- fig %>% add_annotations(x         = 1,
                                 y         = 1,
                                 yref      = 'paper',
                                 yanchor   = "top",
                                 xanchor   = "right",
                                 xref      = "paper",
                                 text      = "<b>Unrelated</b>",
                                 font      = list(size=20),
                                 showarrow = FALSE
                                )
  
  
  
  # Add labels for each individual (Format labels as "ind1 - ind2")
  pretty_labels = gsub('^([MT]{2,3}[0-9]{1,3})([MT]{2,3}[0-9]{1,3})$', '\\1 - \\2', main_only$PairIndividuals)
  fig <- fig %>% add_annotations(x         = main_positions,
                                 y         = main_only$Normalized2AlleleDifference-main_only$StandardError*2,
                                 textangle = 45,
                                 yanchor   = 'top',
                                 xanchor   = 'left',
                                 text      = paste0("<b>",pretty_labels, "</b>"),
                                 showarrow = FALSE,
				 font      = list(size=16)
                                )
  
  # Add annotations and format legend.
  fig <- layout(fig, 
                shapes = rectangles,
                legend = list(orientation = 'h',
                              x           = -0,
                              y           = -0.05,
                              xref        = 'paper',
                              yref        = 'paper',
                              xanchor     = 'left',
                              yanchor     = 'top',
                              font        = list(size=18)
                             ),
                xaxis  = list(showticklabels=F),
                yaxis  = list(title     = list(text="<b> Normalized mean P̄0</b>", font=list(size=18)),
			      font      = list(size=18),
                              titlefont = list(size=15),
                              tickfont  = list(size=16),
                              range     = c(0.5, max(meansP0$Normalized2AlleleDifference+meansP0$StandardError*2)+0.01)
                             ),
                margin = list(l = 5,
                              r = 5,
                              t = 0,
                              b = 0
                             )
               ) %>%
    plotly::config(
      editable    = TRUE,
      displaylogo = FALSE,
      scrollZoom  = TRUE,
      toImageButtonOptions = list(
        format = optargs$format,
        filename = paste0("READ-boxplot")
      )
    )
  
  
  saveWidget(fig, "READ_plot.html", selfcontained=F, libdir="READ_lib_plot")
  fig
}

if (!interactive()) {
  main()
}







