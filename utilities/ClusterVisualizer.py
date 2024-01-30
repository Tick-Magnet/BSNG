import sys, getopt
import seaborn as sns
import pandas as panda
import matplotlib.pyplot as plot


usageString ="""Example: python ClusterVisualizer.py --2D -i inputfile.csv \n
-i filename : input file name \n
-h, --help : print usage \n
--2D : 2D graph \n
--3D : 3D graph \n
-1 dim1 : indicate what dimension should be used for x (high dimensional data) \n
-2 dim2 : indicate what dimension should be used for y (high dimensional data) \n
-3 dim3 : indicate what dimension should be used for z (high dimensional data) \n"""

def usage():
	print(usageString)
	exit(0)

def main():
	inputFile = ""
	twoDim = True
	xDim = 1
	yDim = 2
	zDim = 3
	
	opts, args = getopt.getopt(sys.argv[1:], "hi:1:2:3:", ["3D", "2D", "help"])
	for o, a in opts:
		if o == "-i":
			inputFile = a
		elif o in ("--3D"):
			print("3d")
			twoDim = False
		elif o in ("--2D"):
			print("2D")
			twoDim = True
		elif o in ("--help", "-h"):
			usage()
		elif o == "-1":
			xDim = a
		elif o == "-2":
			yDim = a
		elif o == "-3":
			zDim = a
		else:
			usage()
			exit(1)
		
	
	
	sns.set_theme()
	sns.color_palette("husl");

	#columnNames = ["cluster"];
	#dataset = panda.read_csv(inputFile, names=columnNames)
	dataset = panda.read_csv(inputFile)

	#Changing cluster column to category type so the colors vary more
	dataset.iloc[:,0] = dataset.iloc[:,0].astype('category')

	figure = plot.figure(figsize=(220,220))
	
	if twoDim == True:
		#sns.scatterplot(data=dataset, x= "x", y="y", hue="cluster")
		sns.scatterplot(data=dataset, x= dataset.iloc[:,xDim], y=dataset.iloc[:,yDim], hue=dataset.iloc[:,0])
	else:
		print("3D")
		
		ax = figure.add_subplot(111, projection="3d")
		ax.scatter(xs= dataset.iloc[:,xDim], ys=dataset.iloc[:,yDim], zs = dataset.iloc[:,zDim], c=dataset.iloc[:,0], marker='o');
		ax.set_xlabel("X Axis")
		ax.set_ylabel("Y Axis")
		ax.set_zlabel("Z Axis")
		ax.xaxis.gridlines.set_lw(5.0)
	plot.show()
	
if __name__ == "__main__":
    main()
