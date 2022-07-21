# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006351","transcription, DNA-templated",2.18702440674145,9.80984983294573,0.656408231428349,0,"transcription, DNA-templated"),
                     c("GO:0043170","macromolecule metabolic process",32.1289668000676,10.4281833313121,0.766147755106947,0.20709796,"transcription, DNA-templated"),
                     c("GO:0009058","biosynthetic process",13.618172594555,9.18048327486838,0.823279986117659,0.29181856,"transcription, DNA-templated"),
                     c("GO:0046483","heterocycle metabolic process",14.8526013189786,4.33263526202276,0.80007621273407,0.30014126,"transcription, DNA-templated"),
                     c("GO:0006725","cellular aromatic compound metabolic process",15.1400710219266,4.03873623800409,0.799171636971088,0.30204442,"transcription, DNA-templated"),
                     c("GO:0034641","cellular nitrogen compound metabolic process",18.0034947297221,4.80645938745556,0.780874899228254,0.32040235,"transcription, DNA-templated"),
                     c("GO:1901360","organic cyclic compound metabolic process",16.4985062848768,3.68912229124452,0.803981035735235,0.34957127,"transcription, DNA-templated"),
                     c("GO:0044260","cellular macromolecule metabolic process",25.3311538244744,8.68183813680843,0.740631902214814,0.40932058,"transcription, DNA-templated"),
                     c("GO:0010467","gene expression",10.5067358097063,8.63440148389479,0.798100539340752,0.42501633,"transcription, DNA-templated"),
                     c("GO:0044237","cellular metabolic process",39.9244687447156,9.64977821467881,0.762520967597296,0.44466065,"transcription, DNA-templated"),
                     c("GO:0006807","nitrogen compound metabolic process",36.3113691449186,8.15146609588001,0.770422156122056,0.47283683,"transcription, DNA-templated"),
                     c("GO:0044238","primary metabolic process",39.2593427653458,8.83900693147974,0.764917127621245,0.49275358,"transcription, DNA-templated"),
                     c("GO:0008150","biological_process",100,12.45154581657,1,0,"biological_process"),
                     c("GO:0008152","metabolic process",44.1519643763035,9.32062690148257,1,0,"metabolic process"),
                     c("GO:0009889","regulation of biosynthetic process",23.5161490333127,8.94649885164399,0.767105963184895,0,"regulation of biosynthetic process"),
                     c("GO:0050789","regulation of biological process",65.5205456287695,7.27590047842305,0.903507906755462,0.11857156,"regulation of biosynthetic process"),
                     c("GO:0048519","negative regulation of biological process",29.1471732145877,2.21750011696077,0.912948775168882,0.13394576,"regulation of biosynthetic process"),
                     c("GO:0048522","positive regulation of cellular process",31.5934840200665,1.57099969482502,0.902634435749933,0.14079967,"regulation of biosynthetic process"),
                     c("GO:0019222","regulation of metabolic process",37.489431260921,7.49914241814262,0.903956252813536,0.15796021,"regulation of biosynthetic process"),
                     c("GO:0050794","regulation of cellular process",62.7924017811848,5.1993688432384,0.878973180378872,0.24969315,"regulation of biosynthetic process"),
                     c("GO:0031324","negative regulation of cellular metabolic process",14.0409221577138,2.69173237714189,0.756735980917638,0.35963427,"regulation of biosynthetic process"),
                     c("GO:0010468","regulation of gene expression",26.9939687728989,8.63638713240704,0.74653828537545,0.44489616,"regulation of biosynthetic process"),
                     c("GO:0006355","regulation of transcription, DNA-templated",19.1251902373034,8.79697889623463,0.703086234406766,0.46799282,"regulation of biosynthetic process"),
                     c("GO:0009987","cellular process",85.3277718279691,10.4805799962645,1,0,"cellular process"),
                     c("GO:0032502","developmental process",31.277831012908,2.07981662841049,1,0,"developmental process"),
                     c("GO:0048705","skeletal system morphogenesis",1.26261202863424,1.84578919621559,0.970878711854558,0,"skeletal system morphogenesis"),
                     c("GO:0048731","system development",23.4766924074178,3.0249506812842,0.95178570879929,0.35302259,"skeletal system morphogenesis"),
                     c("GO:0065007","biological regulation",69.5113015049884,6.08999781586737,1,0,"biological regulation"),
                     c("GO:0016043","cellular component organization",29.6150160644834,2.21005002318407,0.992235274543389,0.00923916,"cellular component organization"),
                     c("GO:0071840","cellular component organization or biogenesis",30.7367115720647,1.78644964105139,0.992155305030226,0.01942789,"cellular component organization or biogenesis"),
                     c("GO:0009966","regulation of signal transduction",16.2561298686658,1.80258421093678,0.922217822062209,0.0893218,"regulation of signal transduction"),
                     c("GO:0023051","regulation of signaling",18.3924243278282,1.60990340154366,0.925697690873376,0.09283325,"regulation of signaling"),
                     c("GO:0010646","regulation of cell communication",18.3135110760386,1.66139453487306,0.919304821473897,0.09853671,"regulation of cell communication"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = TRUE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

# dev.off()

