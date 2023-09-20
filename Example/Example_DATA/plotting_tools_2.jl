cs1 = cgrad([:black, :white , :blue]) 
using PyCall
pyimport("numpy")
skl=pyimport("sklearn.decomposition")
PCA=skl.PCA
pyimport("matplotlib")
using PyPlot
plt=PyPlot
grd=pyimport("matplotlib.gridspec")
pca=PCA(n_components=2)

function plot_stat_check(nat_matrix, art_matrix, cm_matrix)
	nat_matrix=reweighted_sample(nat_matrix,12000,0.8)   
	one_hot_encoded_nat=one_hot_encode(nat_matrix,5)
	one_hot_encoded_art=one_hot_encode(art_matrix,5)
	one_hot_encoded_cm=one_hot_encode(cm_matrix,5)
	pca.fit(one_hot_encoded_nat)
	projection_nat=pca.transform(one_hot_encoded_nat)
	projection_art=pca.transform(one_hot_encoded_art)
	projection_cm=pca.transform(one_hot_encoded_cm)
	
	limx1=minimum(projection_nat[:,1])-2
	limx2=maximum(projection_nat[:,1])+2
	
	limy1=minimum(projection_nat[:,2])-2
	limy2=maximum(projection_nat[:,2])+2
	
	fig = plt.figure(figsize=(15, 3))
	gs = grd.GridSpec(nrows=1, ncols=4)

	ax1=fig.add_subplot(get(gs,(0,1)))
	
	plt.hist2D(projection_art[1:end,1],projection_art[1:end,2],bins=30,cmin=1)
	PyPlot.xlim(limx1,limx2)
	PyPlot.ylim(limy1,limy2)
	PyPlot.xlabel("Principal Component 1")
	PyPlot.ylabel("Principal Component 2")
	PyPlot.title("E.A.A. model")

	ax2 = fig.add_subplot(get(gs,(0,0)))

	plt.hist2D(projection_nat[1:end,1],projection_nat[1:end,2],bins=30,cmin=1)
	PyPlot.xlim(limx1,limx2)
	PyPlot.ylim(limy1,limy2)
	PyPlot.xlabel("Principal Component 1")
	PyPlot.ylabel("Principal Component 2")
	PyPlot.title("Natural")


	ax0 = fig.add_subplot(get(gs,(0,-1)))
	
	c1,c2=correlation_comparison_plot_tool(nat_matrix,art_matrix,25000,5)
	plt.scatter(c1,c2)
	scoring=round(cor(correlation_two_point(nat_matrix,5,0),correlation_two_point(art_matrix,5,0));digits=2)
	PyPlot.xlabel("Cij E.E.A.")
	PyPlot.ylabel("Cij Natural")

	PyPlot.title("E.A.A Cij Pearson = $(scoring)")


	ax0 = fig.add_subplot(get(gs,(0,2)))

	plt.hist2D(projection_cm[1:end,1],projection_cm[1:end,2],bins=30,cmin=1)
	PyPlot.xlim(limx1,limx2)
	PyPlot.ylim(limy1,limy2)
	PyPlot.xlabel("Principal Component 1")
	PyPlot.ylabel("Principal Component 2")
	PyPlot.title("Covariance Model")

	ax3 = fig.add_axes([0.92, 0.06, 0.07, 0.35])
	c3,c4=correlation_comparison_plot_tool(nat_matrix,cm_matrix,25000,5)
	scoring=round(cor(correlation_two_point(nat_matrix,5,0),correlation_two_point(cm_matrix,5,0));digits=2)
	plt.scatter(c3,c4, color="red",s=1.5)
	PyPlot.title("CM $scoring")

	PyPlot.subplots_adjust(0,0,1,1,0.5)

	plt.show()
return 
end

