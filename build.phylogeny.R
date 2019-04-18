library(ape)
library(phytools)
#png('figRaxml10ktree.png',width=2000,height=2000,res=400,pointsize=5)
par(xpd=T)							# allow text in the margins

tree<-read.tree("RAxML_bipartitions.UFMG.CM.Y641.mfa")		# read in tree
plot.phylo(tree,cex=0.5, label.offset=0.0003)

#midpoint root the tree until we add in Taiwan
#mtree<-midpoint.root(tree = tree)
#plot.phylo(mtree,cex=0.5, label.offset=0.0002)

rtree<-root(tree, outgroup = c("GE14S01.7B", "EN14S01", "EM14S01.3B"), resolve.root = TRUE)
plot.phylo(rtree,cex=0.5, label.offset=0.0002)
tiplabels()

# rotate the tree so root is at the base									
x<-c(rtree$tip.label[c(12:15)])
rrtree<-rotateConstr(rtree,x)
plot(rrtree, cex =0.4, label.offset=0.0003)

bs<-as.numeric(rrtree$node)				# show bootstraps >= 70%
#bs[bs<70]<-NA

nodelabels(bs,frame="none",cex=0.4,bg="white",adj = c(1.3,-0.2))
add.scale.bar(0.007,18,length = 0.00145,cex=0.8)

# color the clades by population
# function from http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

#set population colors
winecolor<-"#762a83"
paocolor<-"#005a32"
malcolor<-"#308CDF"
sakecolor<-"#e41a1c"
wacolor<-"orange"
taiwancolor<-'#87CEFA'
eucolor<-"#a6dba0"		
ncocolor<-"#5aae61"
nacolor<-"brown"
jap1color<-"#5ddb3e"
jap2acolor<-"#acdb3e"
jap2bcolor<-"#D5EEA8"
jap2ccolor<-"#20aa8d" 
brazil1color<-"#00ac77"
brazil2color<-"#49740c"
brazil3color<-"#5dcc52"
brazil3acolor<-"#577f1c"
#brazil4color<-"#abfa1c"
asianillands1color<-"#84124e"
asianillands2color<-"#db49ee"
asianillands3color<-"#dba8cd"
palmwine<-"#beaa0e"
africanbeer<-"#c95d01"
eastasiacolor<-"gray"
brazilbioethanol<-"tan"
sasoilcolor<-"#817b6a"
beercolor<-"#e67c2b"
philippinescolor<-"#b6cfeb"
njoakcolor<-"#abfa1c"
jap2dcolor<-"#a2ae8c"
bahamascolor<-"#7ac0c0"
philfruitcolor<-"#d64384"
philfruit2color<-"#893d60"

#set ecological category colors
fermentcolor<-"black"
treecolor<-"#12b216"
fruit.flowercolor<-"#6cccd0"			
soil.unknowncolor<-"#865b0d"

#set population nodes 
nodelabels()

winenode<-60
paonode<-110
malnode<-81
sakenode<-90
wanode<-77
tawainnode<-74
eunode<-64
nconode<-107
jap1node<-112
jap2anode<-104
jap2bnode<-100
brazil1node<-84
brazil3node<-70
philfruitnode<-95
  
##note to self the order that this is read into R will determine if all colors are seen. 
pew<-2 #set population edge width
ec<-rep("black",length(rrtree$edge[,2]))	# default edge color = black
ew<-rep(2,length(rrtree$edge[,2]))		# default edge width = 1 makes black and color lines same size
tc<-rep("black",length(rrtree$tip))		# default tip color = black 

#tawain population
d<-getDescendants(rrtree,node=tawainnode)		
for (x in d) { ec[rrtree$edge[,2]==x]<-taiwancolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-taiwancolor}
ew[rrtree$edge[,2]==tawainnode]<-pew

# wa population
d<-getDescendants(rrtree,node=wanode)	
for (x in d) { ec[rrtree$edge[,2]==x]<-wacolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-wacolor }
ec[rrtree$edge[,2]==wanode]<-wacolor	# color the edge leading to this clade
ew[rrtree$edge[,2]==wanode]<-pew

# euo population
d<-getDescendants(rrtree,node=eunode)		
for (x in d) { ec[rrtree$edge[,2]==x]<-eucolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-eucolor }
ec[rrtree$edge[,2]==eunode]<-eucolor
ew[rrtree$edge[,2]==eunode]<-pew

# wine population
d<-getDescendants(rrtree,node=winenode)		
for (x in d) { ec[rrtree$edge[,2]==x]<-winecolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-winecolor }
ec[rrtree$edge[,2]==winenode]<-winecolor		# color the edge leading to this clade

#brazil 3 population
d<-getDescendants(rrtree,node=brazil3node)	
for (x in d) { ec[rrtree$edge[,2]==x]<-brazil3color; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-brazil3color }
ec[rrtree$edge[,2]==brazil3node]<-brazil3color
ew[rrtree$edge[,2]==brazil3node]<-pew

#brazil 1 population
d<-getDescendants(rrtree,node=brazil1node)	
for (x in d) { ec[rrtree$edge[,2]==x]<-brazil1color; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-brazil1color }
ec[rrtree$edge[,2]==brazil1node]<-brazil1color
ew[rrtree$edge[,2]==brazil1node]<-pew

# mal population
d<-getDescendants(rrtree,node=malnode)	
for (x in d) { ec[rrtree$edge[,2]==x]<-malcolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-malcolor }
ec[rrtree$edge[,2]==malnode]<-malcolor	# color the edge leading to this clade
ew[rrtree$edge[,2]==malnode]<-pew

# sake population
d<-getDescendants(rrtree,node=sakenode)	
for (x in d) { ec[rrtree$edge[,2]==x]<-sakecolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-sakecolor }
ec[rrtree$edge[,2]==sakenode]<-sakecolor
ew[rrtree$edge[,2]==sakenode]<-pew

#philippines fruit population
d<-getDescendants(rrtree,node=philfruitnode)		
for (x in d) { ec[rrtree$edge[,2]==x]<-philfruitcolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-philfruitcolor}
ew[rrtree$edge[,2]==philfruitnode]<-pew
ew[rrtree$edge[,2]==philfruitnode]<-pew

# japan population 2a
d<-getDescendants(rrtree,node=jap2anode)	
for (x in d) { ec[rrtree$edge[,2]==x]<-jap2acolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-jap2acolor }
ec[rrtree$edge[,2]==jap2anode]<-jap2acolor

# japan population 2b
d<-getDescendants(rrtree,node=jap2bnode)	
for (x in d) { ec[rrtree$edge[,2]==x]<-jap2bcolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-jap2bcolor }
ec[rrtree$edge[,2]==jap2bnode]<-jap2bcolor

# nco population
d<-getDescendants(rrtree,node=nconode)	
for (x in d) { ec[rrtree$edge[,2]==x]<-ncocolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-ncocolor }
ec[rrtree$edge[,2]==nconode]<-ncocolor
ew[rrtree$edge[,2]==nconode]<-pew

# japan population 1
d<-getDescendants(rrtree,node=jap1node)	
for (x in d) { ec[rrtree$edge[,2]==x]<-jap1color; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-jap1color }
ec[rrtree$edge[,2]==jap1node]<-jap1color

# pao population
d<-getDescendants(rrtree,node=paonode)	
for (x in d) { ec[rrtree$edge[,2]==x]<-paocolor; ew[rrtree$edge[,2]==x]<-pew; tc[x]<-paocolor }
ec[rrtree$edge[,2]==paonode]<-paocolor

tc[grep("UFMG.CM.Y639",rrtree$tip)]<-brazil2color
tc[grep("UFMG.CM.Y642",rrtree$tip)]<-brazil3color
tc[grep("UFMG.CM.Y263",rrtree$tip)]<-brazil3acolor
#tc[grep("UFMG.CM.Y260",rrtree$tip)]<-brazil4color
tc[grep("YJM1463",rrtree$tip)]<-nacolor
tc[grep("ZP793",rrtree$tip)]<-jap2ccolor
tc[grep("S8BM.32.4D.a",rrtree$tip)]<-asianillands1color
tc[grep("CBS1594",rrtree$tip)]<-asianillands2color
tc[grep("CBS3000",rrtree$tip)]<-asianillands3color
tc[grep("CBS1419",rrtree$tip)]<-palmwine
tc[grep("CBS4456",rrtree$tip)]<-africanbeer
tc[grep("CLIB219.2b",rrtree$tip)]<-eastasiacolor
tc[grep("CBS7959",rrtree$tip)]<-brazilbioethanol
tc[grep("CBS2888.1b",rrtree$tip)]<-sasoilcolor
tc[grep("CBS1463",rrtree$tip)]<-beercolor
tc[grep("Y10.1b",rrtree$tip)]<-philippinescolor
tc[grep("YPS1000",rrtree$tip)]<-njoakcolor
tc[grep("ZP780",rrtree$tip)]<-jap2dcolor
tc[grep("UWOPS83.787.3",rrtree$tip)]<-bahamascolor
tc[grep("YJM1401",rrtree$tip)]<-philfruit2color
tc[grep("GE14S01.7B",rrtree$tip)]<-taiwancolor

plot.phylo(rrtree,cex=.9,edge.color=ec,tip.color=tc,edge.width=ew,font=2, label.offset=0.0001)

nodelabels(bs,frame="none",cex=0.8,bg="white",adj = c(1.2,-0.2))
add.scale.bar(0,-2, length = 0.00145, cex=0.8)

##add dots to represent ecological categoreis 
#trees
tiplabels("",grep("ZP",rrtree$tip),pch=20,frame="none",col=treecolor,adj=0.5)
tiplabels("",grep("YPS",rrtree$tip),pch=20,frame="none",col=treecolor,adj=0.5)
tiplabels("",grep("SD",rrtree$tip),pch=20,frame="none",col=treecolor,adj=0.5)
tiplabels("",grep("SM",rrtree$tip),pch=20,frame="none",col=treecolor,adj=0.5)
tiplabels("",grep("S8BM.32.4D.a	",rrtree$tip),pch=20,frame="none",col=treecolor,adj=0.5)
tiplabels("",grep("YJM1273",rrtree$tip),pch=20,frame="none",col=treecolor,adj=0.5)

#Fruit or Flower
tiplabels("",grep("UWOPS",rrtree$tip),pch=20,frame="none",col=fruit.flowercolor,adj=0.5)
tiplabels("",grep("UFMG",rrtree$tip),pch=20,frame="none",col=fruit.flowercolor,adj=0.5)
tiplabels("",grep("Y10",rrtree$tip),pch=20,frame="none",col=fruit.flowercolor,adj=0.5)
tiplabels("",grep("YJM140",rrtree$tip),pch=20,frame="none",col=fruit.flowercolor,adj=0.5)
tiplabels("",grep("YJM1479",rrtree$tip),pch=20,frame="none",col=fruit.flowercolor,adj=0.5)

#Fermentation
tiplabels("",grep("Y12",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("K11",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("YJM1388",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("DBVPG6044",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("YJM1439",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("YJM1248",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("L1374",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("BC187",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("YJM1332",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("CBS3000",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("CBS1419",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("CBS4456",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("CLIB219.2b",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)
tiplabels("",grep("CBS1463",rrtree$tip),pch=20,frame="none",col=fermentcolor,adj=0.5)

#Soil and Other
tiplabels("",grep("EN14S01",rrtree$tip),pch=20,frame="none",col=soil.unknowncolor,adj=0.5)
tiplabels("",grep("EM14S01.3B",rrtree$tip),pch=20,frame="none",col=soil.unknowncolor,adj=0.5)
tiplabels("",grep("GE14S01.7B",rrtree$tip),pch=20,frame="none",col=soil.unknowncolor,adj=0.5)
tiplabels("",grep("YJM1463",rrtree$tip),pch=20,frame="none",col=soil.unknowncolor,adj=0.5)
tiplabels("",grep("CBS1594",rrtree$tip),pch=20,frame="none",col=soil.unknowncolor,adj=0.5)
tiplabels("",grep("CBS7959",rrtree$tip),pch=20,frame="none",col=soil.unknowncolor,adj=0.5)
tiplabels("",grep("CBS2888.1b",rrtree$tip),pch=20,frame="none",col=soil.unknowncolor,adj=0.5)


#legend.coord<-locator(1)
#legend.coord

lineagenames<-c( "Brazil 3a Oak", "Beer", "Wine", "European Oak", "African Beer", "North Africa",
                "Brazilian Bioethanol", "Brazil 3 Oak",  "Brazil 2 Oak", "Malaysia", "New Jersey Oak", 
                "Brazil 1 Oak", "Japan 2d Oak", "Japan 2a Oak", "Japan 2c Oak", "Japan 2b Oak", "Pennsylvania Oak",
                "Japan 1 Oak", "North Carolina Oak",  "Bahamas",  "Sake", "Philippines", "Asian Islands 2",
                "Asian Islands 1", "Philippines Fruit 1", "Philippines Fruit 2",  "Asian Islands 3", 
                "South African Soil",  "West AFrica",  "African Palm Wine", "Far East Asia",
                "Taiwan")

lineagecolors<-c( brazil3acolor, beercolor, winecolor, eucolor, africanbeer,  nacolor,
                 brazilbioethanol, brazil3color, brazil2color, malcolor, njoakcolor,
                 brazil1color, jap2dcolor, jap2acolor, jap2ccolor, jap2bcolor, paocolor,
                 jap1color, ncocolor, bahamascolor, sakecolor, philippinescolor, asianillands2color,
                 asianillands1color, philfruitcolor, philfruit2color, asianillands3color,  
                 sasoilcolor,  wacolor, palmwine,eastasiacolor,
                 taiwancolor)

lineagelty<-rep(1,7)
lineagelwd<-rep(2,7)
legend(0.001,60.5, lineagenames,col=lineagecolors,lty=lineagelty,lwd=lineagelwd,bty="n",cex=1,title="Reference Populations",title.adj=0,text.col=lineagecolors,title.col="black", text.font = 2)

ecologicalnames<-c("Fermentations", "Fruit or Flower", "Tree", "Soil or Other")
ecologicalcolors<-c(fermentcolor, fruit.flowercolor, treecolor, soil.unknowncolor)
legend(0.005, 60.5,ecologicalnames,col=ecologicalcolors,bty="n",cex=1,title = "Ecological Caterogies",title.adj=0,text.col=ecologicalcolors,title.col="black",text.font = 2, pch =19)
