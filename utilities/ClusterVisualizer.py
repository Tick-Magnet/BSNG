import sys
import seaborn as sns
import pandas as panda
import matplotlib.pyplot as plot
from matplotlib.widgets import Slider

sns.set_theme()


columnNames = ["cluster", "x", "y"];
dataset = panda.read_csv(sys.argv[1], names=columnNames)

figure = plot.figure(figsize=(220,220))


sns.scatterplot(data=dataset, x= "x", y="y", hue="cluster")



plot.show()
