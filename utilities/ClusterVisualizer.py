import seaborn as sns
import pandas as panda
import matplotlib.pyplot as plot
from matplotlib.widgets import Slider

sns.set_theme()

#See https://www.tutorialspoint.com/scroll-backwards-and-forwards-through-matplotlib-plots for scrolling method

columnNames = ["cluster", "x", "y"];
dataset = panda.read_csv("test.csv", names=columnNames)

figure = plot.figure(figsize=(20,20))


sns.scatterplot(data=dataset, x= "x", y="y", hue="cluster")



plot.show()
