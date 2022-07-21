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
revigo.data <- rbind(c("GO:0001738","morphogenesis of a polarized epithelium",0.282661541070722,1.78433236782293,0.999135063930978,0,"morphogenesis of a polarized epithelium"),
                     c("GO:0002478","antigen processing and presentation of exogenous peptide antigen",0.214822771213749,4.4824366634185,0.968905070262002,0,"antigen processing and presentation of exogenous peptide antigen"),
                     c("GO:0060218","hematopoietic stem cell differentiation",0.090451693142631,1.80758924045032,0.981570371904805,0.26655813,"antigen processing and presentation of exogenous peptide antigen"),
                     c("GO:0019882","antigen processing and presentation",0.531403697212957,2.28511964709216,0.985658005370229,0.30803218,"antigen processing and presentation of exogenous peptide antigen"),
                     c("GO:0008150","biological_process",100,20.5433505008239,1,0,"biological_process"),
                     c("GO:0008152","metabolic process",43.8181920967833,35.8009227930301,1,0,"metabolic process"),
                     c("GO:0009987","cellular process",85.1263497088586,21.0835936349364,1,0,"cellular process"),
                     c("GO:0010608","posttranscriptional regulation of gene expression",3.0131720278139,14.265396693887,0.918165238777769,0,"posttranscriptional regulation of gene expression"),
                     c("GO:0006521","regulation of cellular amino acid metabolic process",0.0734920006783877,4.50733168327864,0.948897273530738,0.11574912,"posttranscriptional regulation of gene expression"),
                     c("GO:0033238","regulation of cellular amine metabolic process",0.169596924642433,2.26262418993652,0.944477318412534,0.12554452,"posttranscriptional regulation of gene expression"),
                     c("GO:1903311","regulation of mRNA metabolic process",1.62247724574594,6.84521134727866,0.920397797032096,0.18352065,"posttranscriptional regulation of gene expression"),
                     c("GO:0062125","regulation of mitochondrial gene expression",0.158290462999604,2.34743751479664,0.942897162222546,0.19325517,"posttranscriptional regulation of gene expression"),
                     c("GO:0045815","positive regulation of gene expression, epigenetic",0.175250155463848,2.46144632664236,0.933946393356934,0.19525236,"posttranscriptional regulation of gene expression"),
                     c("GO:0031329","regulation of cellular catabolic process",4.57911696534569,5.84492837361737,0.908700165344318,0.20187102,"posttranscriptional regulation of gene expression"),
                     c("GO:0009894","regulation of catabolic process",5.3649160495223,7.37890992774562,0.922222369306163,0.206457,"posttranscriptional regulation of gene expression"),
                     c("GO:0060968","regulation of gene silencing",0.390072926677596,2.02611088328014,0.935987469971949,0.21251716,"posttranscriptional regulation of gene expression"),
                     c("GO:0034248","regulation of cellular amide metabolic process",2.74747017920742,4.94937950137051,0.92157606168067,0.21779321,"posttranscriptional regulation of gene expression"),
                     c("GO:0061418","regulation of transcription from RNA polymerase II promoter in response to hypoxia",0.0734920006783877,6.79111960801894,0.898036299035755,0.21956837,"posttranscriptional regulation of gene expression"),
                     c("GO:0040029","regulation of gene expression, epigenetic",0.729266775962463,5.38425837305549,0.93237339234843,0.22830396,"posttranscriptional regulation of gene expression"),
                     c("GO:0010033","response to organic substance",15.0149810616767,2.76761537495874,0.98170856118114,0.2514283,"posttranscriptional regulation of gene expression"),
                     c("GO:0006417","regulation of translation",2.43088925320821,5.75967996881783,0.8911409897049,0.26637964,"posttranscriptional regulation of gene expression"),
                     c("GO:0043484","regulation of RNA splicing",0.836678161569337,1.57019880862467,0.923076068190404,0.27954208,"posttranscriptional regulation of gene expression"),
                     c("GO:0051246","regulation of protein metabolic process",14.2178755158573,7.40585559537935,0.896132836919644,0.28737984,"posttranscriptional regulation of gene expression"),
                     c("GO:0010629","negative regulation of gene expression",4.90135112216632,14.2271643621483,0.877768011661335,0.29503942,"posttranscriptional regulation of gene expression"),
                     c("GO:0006283","transcription-coupled nucleotide-excision repair",0.0565323082141444,1.90469903066537,0.828020260172992,0.30726982,"posttranscriptional regulation of gene expression"),
                     c("GO:0045910","negative regulation of DNA recombination",0.25439538696365,1.30611908793955,0.911890123687149,0.35426348,"posttranscriptional regulation of gene expression"),
                     c("GO:0009892","negative regulation of metabolic process",16.4452484594946,10.2220231713853,0.880601715889951,0.36582086,"posttranscriptional regulation of gene expression"),
                     c("GO:0016458","gene silencing",0.831024930747922,2.72416728895381,0.893157041021292,0.40798761,"posttranscriptional regulation of gene expression"),
                     c("GO:0031396","regulation of protein ubiquitination",1.16456554921137,1.77178468499161,0.919602334995572,0.42779096,"posttranscriptional regulation of gene expression"),
                     c("GO:0060255","regulation of macromolecule metabolic process",34.4168692407711,8.40475955149823,0.884371434960132,0.43199703,"posttranscriptional regulation of gene expression"),
                     c("GO:0006974","cellular response to DNA damage stimulus",4.03075357566849,6.46662208183846,0.957382394313485,0.43300999,"posttranscriptional regulation of gene expression"),
                     c("GO:0042769","DNA damage response, detection of DNA damage",0.0734920006783877,1.34747902357484,0.968559975507352,0.43300999,"posttranscriptional regulation of gene expression"),
                     c("GO:0070482","response to oxygen levels",1.84860647860252,3.09230709265733,0.978773659169685,0.46363477,"posttranscriptional regulation of gene expression"),
                     c("GO:0033554","cellular response to stress",8.48549946294307,15.7544618493379,0.962204724475793,0.47963796,"posttranscriptional regulation of gene expression"),
                     c("GO:2000112","regulation of cellular macromolecule biosynthetic process",22.1889309740517,2.79096080119664,0.873304633926957,0.48230717,"posttranscriptional regulation of gene expression"),
                     c("GO:0009889","regulation of biosynthetic process",23.5683192944768,2.65705850737022,0.894983893351842,0.49388724,"posttranscriptional regulation of gene expression"),
                     c("GO:0016032","viral process",1.43026739781785,24.3614819802361,1,0,"viral process"),
                     c("GO:0019080","viral gene expression",0.277008310249307,22.5376372367054,0.996324446128815,0,"viral gene expression"),
                     c("GO:0034622","cellular protein-containing complex assembly",4.42082650234609,30.0416864901722,0.892395166960955,0,"cellular protein-containing complex assembly"),
                     c("GO:0022411","cellular component disassembly",1.57159816835321,3.37145388337854,0.941897076748827,0.23465096,"cellular protein-containing complex assembly"),
                     c("GO:0007005","mitochondrion organization",2.30086494431568,13.0046947658058,0.924542785007159,0.24763174,"cellular protein-containing complex assembly"),
                     c("GO:0061024","membrane organization",4.37560065577478,2.26015808275786,0.933829998556759,0.27310681,"cellular protein-containing complex assembly"),
                     c("GO:0032200","telomere organization",0.45791169653457,4.06113621073378,0.922754226707243,0.29582015,"cellular protein-containing complex assembly"),
                     c("GO:0048285","organelle fission",1.87121940188818,2.77293930487855,0.926388096842109,0.34954815,"cellular protein-containing complex assembly"),
                     c("GO:0016043","cellular component organization",29.7077279665329,32.7370038823907,0.913563563756314,0.3760284,"cellular protein-containing complex assembly"),
                     c("GO:0051276","chromosome organization",5.7549889761999,10.1089542019362,0.915102442415584,0.40880723,"cellular protein-containing complex assembly"),
                     c("GO:0043933","protein-containing complex subunit organization",7.54140991576686,25.7721585065052,0.928538585192123,0.42889898,"cellular protein-containing complex assembly"),
                     c("GO:0042254","ribosome biogenesis",1.63943693821019,18.351970386497,0.909033223922264,0.46579103,"cellular protein-containing complex assembly"),
                     c("GO:0022613","ribonucleoprotein complex biogenesis",2.38001017581548,21.6408744764468,0.920222824567881,0.49111213,"cellular protein-containing complex assembly"),
                     c("GO:0051179","localization",28.7127593419639,2.57412417401352,1,0,"localization"),
                     c("GO:0006613","cotranslational protein targeting to membrane",0.118717847249703,16.7772526649581,0.913955887862055,0.00486079,"cotranslational protein targeting to membrane"),
                     c("GO:0006839","mitochondrial transport",1.15891231838996,4.73148685948869,0.948728319045326,0.21765143,"cotranslational protein targeting to membrane"),
                     c("GO:0051641","cellular localization",12.7254225790039,12.4525860575054,0.931184060121942,0.27289693,"cotranslational protein targeting to membrane"),
                     c("GO:0006403","RNA localization",1.03454124031884,1.98211933367975,0.933740458180632,0.31531678,"cotranslational protein targeting to membrane"),
                     c("GO:0071705","nitrogen compound transport",8.85861269715642,7.26691866448147,0.931568087252113,0.39664369,"cotranslational protein targeting to membrane"),
                     c("GO:0033036","macromolecule localization",12.9119791961106,7.95225949909517,0.938797021917041,0.43301671,"cotranslational protein targeting to membrane"),
                     c("GO:0051169","nuclear transport",1.19848493413986,2.0102667488078,0.916143907281034,0.44226858,"cotranslational protein targeting to membrane"),
                     c("GO:0072594","establishment of protein localization to organelle",1.67900955396009,12.587239024165,0.894573578551469,0.45193914,"cotranslational protein targeting to membrane"),
                     c("GO:0051668","localization within membrane",2.83792187235005,4.92211543293267,0.917128648413106,0.4704432,"cotranslational protein targeting to membrane"),
                     c("GO:0006457","protein folding",1.1758720108542,1.80633132335012,0.993586353918246,0.00633472,"protein folding"),
                     c("GO:0061919","process utilizing autophagic mechanism",1.53202555260331,1.65643650029214,0.993374876815514,0.00656439,"process utilizing autophagic mechanism"),
                     c("GO:0007059","chromosome segregation",1.58855786081746,3.61173204106109,0.993344789260699,0.00659715,"chromosome segregation"),
                     c("GO:0006412","translation",2.00689694160213,27.510164146586,0.743971816622598,0.0068166,"translation"),
                     c("GO:0046034","ATP metabolic process",1.1362993951043,4.000341289956,0.916243734516576,0.11055793,"translation"),
                     c("GO:0009060","aerobic respiration",0.864944315676409,8.05936784516778,0.868864929033136,0.12884573,"translation"),
                     c("GO:0006091","generation of precursor metabolites and energy",2.1934535587088,6.19825287594235,0.892052269403836,0.14436015,"translation"),
                     c("GO:0044237","cellular metabolic process",39.6687206738651,44.4562910909582,0.83973899971842,0.19542101,"translation"),
                     c("GO:0006520","cellular amino acid metabolic process",1.54333201424614,3.59386969756245,0.850634339413447,0.23637028,"translation"),
                     c("GO:0044281","small molecule metabolic process",8.72293515744248,4.99192460967878,0.886743600590424,0.27561982,"translation"),
                     c("GO:0006793","phosphorus metabolic process",10.1249364011533,4.55320604463459,0.863233246685986,0.28839829,"translation"),
                     c("GO:0009056","catabolic process",10.6167674826163,17.6328348526358,0.882687894507863,0.29271726,"translation"),
                     c("GO:0009058","biosynthetic process",13.5790604330375,20.4144574116319,0.877149442586623,0.3173766,"translation"),
                     c("GO:0046483","heterocycle metabolic process",14.7379727514274,19.1162667402964,0.853465406909606,0.32653104,"translation"),
                     c("GO:0006725","cellular aromatic compound metabolic process",15.0828198315337,17.7059679899362,0.852814775188306,0.32921276,"translation"),
                     c("GO:1901360","organic cyclic compound metabolic process",16.3943693821019,17.7454555415329,0.861864716341899,0.33925737,"translation"),
                     c("GO:0034655","nucleobase-containing compound catabolic process",1.58290462999604,16.2974344521539,0.764924769585652,0.35056933,"translation"),
                     c("GO:0034641","cellular nitrogen compound metabolic process",17.9094352422409,30.0511184676175,0.832294326413329,0.35059692,"translation"),
                     c("GO:0043603","cellular amide metabolic process",4.31341511673922,22.6915855831392,0.813961578732817,0.40040394,"translation"),
                     c("GO:0070647","protein modification by small protein conjugation or removal",4.78828650573803,2.7648609742611,0.817242294922347,0.40052964,"translation"),
                     c("GO:0031145","anaphase-promoting complex-dependent catabolic process",0.107411385606874,7.83849865932845,0.837472891881735,0.40060606,"translation"),
                     c("GO:0044260","cellular macromolecule metabolic process",25.0946916162587,21.2885504011177,0.803340525388725,0.40185837,"translation"),
                     c("GO:0010467","gene expression",10.3510656340098,15.8085724775899,0.843528821066446,0.4293132,"translation"),
                     c("GO:1901564","organonitrogen compound metabolic process",27.022443326361,22.3588806521924,0.825369497376588,0.43308452,"translation"),
                     c("GO:0140053","mitochondrial gene expression",0.497484312284471,7.71254555663999,0.869817320143462,0.44145094,"translation"),
                     c("GO:0043170","macromolecule metabolic process",31.7315846005992,21.0675905916553,0.840758692193378,0.44735906,"translation"),
                     c("GO:0009303","rRNA transcription",0.118717847249703,4.62594822238045,0.804443918214022,0.45819339,"translation"),
                     c("GO:0006807","nitrogen compound metabolic process",35.9319351009102,31.5207199767676,0.848373013073835,0.47590559,"translation"),
                     c("GO:0043412","macromolecule modification",16.0721352252812,9.55281338335334,0.830415034957703,0.48814861,"translation"),
                     c("GO:0044238","primary metabolic process",38.8263892814744,33.2951120037051,0.845454434137082,0.49561512,"translation"),
                     c("GO:0051301","cell division",2.79269602577873,5.30725267232866,0.992837038294262,0.00715292,"cell division"),
                     c("GO:1903047","mitotic cell cycle process",2.84357510317146,14.1050381140334,0.945097283884182,0.00717225,"mitotic cell cycle process"),
                     c("GO:0071840","cellular component organization or biogenesis",30.8101079767087,42.7409802368363,0.989296184016113,0.01115019,"cellular component organization or biogenesis"),
                     c("GO:0007049","cell cycle",6.87998190966137,15.7718055389001,0.99183391451153,0.01242984,"cell cycle"),
                     c("GO:0031647","regulation of protein stability",1.6846627847815,2.27081297256644,0.978712138391189,0.03000366,"regulation of protein stability"),
                     c("GO:1901532","regulation of hematopoietic progenitor cell differentiation",0.180903386285262,3.86659689788049,0.972860683222461,0.03041435,"regulation of hematopoietic progenitor cell differentiation"),
                     c("GO:0090175","regulation of establishment of planar polarity",0.322234156820623,2.24979424998361,0.971134449119157,0.29093563,"regulation of hematopoietic progenitor cell differentiation"),
                     c("GO:1902036","regulation of hematopoietic stem cell differentiation",0.0791452314998021,3.14853782188146,0.974771542527885,0.33776283,"regulation of hematopoietic progenitor cell differentiation"),
                     c("GO:2000736","regulation of stem cell differentiation",0.327887387642037,2.79254876229192,0.971294307143832,0.37765336,"regulation of hematopoietic progenitor cell differentiation"),
                     c("GO:0038093","Fc receptor signaling pathway",0.248742156142235,4.61421360532595,0.943815532858763,0.03143403,"Fc receptor signaling pathway"),
                     c("GO:0006950","response to stress",18.7235004805246,3.42084292890144,0.983106179602061,0.11368993,"Fc receptor signaling pathway"),
                     c("GO:0097193","intrinsic apoptotic signaling pathway",0.819718469105093,2.09036767633648,0.955628665259167,0.18030448,"Fc receptor signaling pathway"),
                     c("GO:0009628","response to abiotic stimulus",6.25812651930578,1.62766612085893,0.985524902646779,0.19615032,"Fc receptor signaling pathway"),
                     c("GO:0090263","positive regulation of canonical Wnt signaling pathway",0.582282774605687,2.67554664175951,0.948911080165641,0.21870843,"Fc receptor signaling pathway"),
                     c("GO:0048522","positive regulation of cellular process",31.251059980779,5.39238338111869,0.943988629038904,0.22604981,"Fc receptor signaling pathway"),
                     c("GO:0080135","regulation of cellular response to stress",3.56718864831251,1.98758780910557,0.95487822718234,0.28739827,"Fc receptor signaling pathway"),
                     c("GO:0070498","interleukin-1-mediated signaling pathway",0.113064616428289,2.85864987129045,0.958481438240408,0.29988107,"Fc receptor signaling pathway"),
                     c("GO:0060070","canonical Wnt signaling pathway",0.440952004070326,1.30639166284386,0.952434679985423,0.33561779,"Fc receptor signaling pathway"),
                     c("GO:0030111","regulation of Wnt signaling pathway",1.79207417038838,1.46854892233098,0.958240842566738,0.33636761,"Fc receptor signaling pathway"),
                     c("GO:0009966","regulation of signal transduction",15.9647238396744,1.43523174360297,0.943118187914158,0.47862136,"Fc receptor signaling pathway"),
                     c("GO:0051338","regulation of transferase activity",4.90135112216632,2.5469334064591,0.975372357914425,0.03491975,"regulation of transferase activity"),
                     c("GO:0010564","regulation of cell cycle process",3.48804341681271,9.42830416363332,0.925214482236198,0.04353772,"regulation of cell cycle process"),
                     c("GO:0033044","regulation of chromosome organization",1.57725139917463,3.48252003571517,0.971312735535961,0.0467518,"regulation of chromosome organization"),
                     c("GO:0019222","regulation of metabolic process",37.2265249590141,9.7358039934942,0.957079688925818,0.06873468,"regulation of metabolic process"),
                     c("GO:0010646","regulation of cell communication",18.073378936062,2.26860247010631,0.960483434232298,0.11063185,"regulation of metabolic process"),
                     c("GO:0023051","regulation of signaling",18.1525241675618,1.9786832437451,0.964139020194444,0.11081125,"regulation of metabolic process"),
                     c("GO:0048519","negative regulation of biological process",28.9219288823563,6.1964488957902,0.959880882378729,0.13397022,"regulation of metabolic process"),
                     c("GO:0048518","positive regulation of biological process",33.9646107750579,4.70026310357757,0.958144187867658,0.144382,"regulation of metabolic process"),
                     c("GO:0051726","regulation of cell cycle",5.60800497484312,7.35014558512478,0.968520514607239,0.07717043,"regulation of cell cycle"),
                     c("GO:0051128","regulation of cellular component organization",13.3472779693595,3.24731880048346,0.962966698688675,0.09946079,"regulation of cellular component organization"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
# pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

