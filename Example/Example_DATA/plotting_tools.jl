using Distances
using StatsBase
using Random
using Statistics
using Plots
using PlotThemes
theme(:mute)
using Colors, ColorSchemes




function correlation_comparison_plot_tool(matrix_1,matrix_2,Number,q)
	correlation_list_1=correlation_two_point(matrix_1,q,0)
	correlation_list_2=correlation_two_point(matrix_2,q,0)
	len=length(correlation_list_1)
	correlation_to_plot_1=zeros(Float32,Number)
	correlation_to_plot_2=zeros(Float32,Number)
	for i in 1:Number
	    j=rand(1:len)
	    correlation_to_plot_1[i]=correlation_list_1[j]
	    correlation_to_plot_2[i]=correlation_list_2[j]
	end
	return correlation_to_plot_1,correlation_to_plot_2
end



function secondary_structure_plot_tools(filename)
	f = open(filename)                            
	lines=readlines(f)
	copy=lines
	file_copy=[]
	for i in 1:length(copy)
	    file_copy=push!(file_copy,[])
	end
	for i in 1:length(copy)
		file_copy[i]=split(copy[i],[' ','[',']',','])
	end
	starting_condition=0
	i=1
	while starting_condition==0
	    if file_copy[i][1]=="HMM"
		starting_condition=1
	    else
	    i=i+1
	    end
	end 
	file_copy=file_copy[i+3:end-1]
	dot_bracket_ss=[]
	molecule=[]
	for i in 1:length(file_copy)
	    if i%3==0
		push!(dot_bracket_ss,file_copy[i][end])
		push!(molecule,file_copy[i][end-3])
	    end
	end
	for i in 1:length(dot_bracket_ss)
	    if dot_bracket_ss[i]=="{" || dot_bracket_ss[i]=="<" || dot_bracket_ss[i]=="[" || dot_bracket_ss[i]=="("
	       
		dot_bracket_ss[i]="("
	    elseif dot_bracket_ss[i]=="}" || dot_bracket_ss[i]==">" || dot_bracket_ss[i]=="]" || dot_bracket_ss[i]==")"
	       
		dot_bracket_ss[i]=")"
	    else
	       
		dot_bracket_ss[i]="."
	    end
	end
	dot_braket_ss=join(dot_bracket_ss)
	len=length(dot_bracket_ss)
	contact_matrix=zeros(Int8,len,len)  
	number_of_contacts=length(dot_bracket_ss[dot_bracket_ss.=="("])
	contact_list=zeros(number_of_contacts,2)                                                 

	ss_proxy=zeros(Int8,len)  

	for i in 1:len
	    if dot_bracket_ss[i]=="("
		ss_proxy[i]=1
	    elseif  dot_bracket_ss[i]==")"
		ss_proxy[i]=2
	    else
		ss_proxy[i]==0
	    end
	end
	contact_matrix=zeros(Int8,len,len)  
	contact_list=zeros(Int64,number_of_contacts,2)  
	idx=1
	for i in 1:len
	    k=0
	    if ss_proxy[i]==2
		ss_proxy[i]=0
		for j in 1:i
		    
		    if k==0
		    if j<i
		        
		        if ss_proxy[i-j]==1
		            contact_matrix[i-j,i]=1
		            
		            contact_list[idx,:]=[i-j,i]
		            idx+=1
		            ss_proxy[i-j]=0
		            k=1
		            
		        end
		    end
		    end
		    
		end
	    end
	end
	return contact_list, contact_matrix ,join(molecule), join(dot_bracket_ss), len
end



function tertiary_plot_tools(len,filename,limit_contact_distance)

f = open(filename)                       
lines=readlines(f) 

file_copy=[]
copy=lines

for i in 1:length(copy)
    file_copy=push!(file_copy,[])
end

for i in 1:length(copy)    
file_copy[i]=split(copy[i],[' ','[',']',','])     
end

distance_matrix=zeros(len,len)
for i in 2:length(file_copy)
    distance_matrix[parse.(Int64,file_copy[i][3]),parse.(Int64,file_copy[i][7])]=parse.(Float64,file_copy[i][10])
end

contact_matrix=zeros(len,len)
contact_list=zeros(Int8,0,2)
for i in 1:len
    for j in i+4:len
        if distance_matrix[i,j]<limit_contact_distance && distance_matrix[i,j]!=0                      
            contact_matrix[i,j]=distance_matrix[i,j]
            contact_list=vcat(contact_list,[i,j]')
        end
    end
end

return contact_list, contact_matrix
end



function edge_interpretation_plot(len,ss_contact_matrix,tertiary_contact_matrix,edge_list)

sec_proxy_list=findall(!iszero, ss_contact_matrix)
proxy_idx_1 = getindex.(sec_proxy_list, 1)
proxy_idx_2 = getindex.(sec_proxy_list, 2)
ss_contact_list=hcat(proxy_idx_1,proxy_idx_2)

ter_proxy_list=findall(!iszero, tertiary_contact_matrix)
proxy_idx_1 = getindex.(ter_proxy_list, 1)
proxy_idx_2 = getindex.(ter_proxy_list, 2)
tertiary_contact_list=hcat(proxy_idx_1,proxy_idx_2)

bar2D=0
bar3D=0
barNeigh=0
barTot=0

proxy_1=[]
proxy_2=[]
proxy_3=[]
proxy_4=[]
sec_edge=zeros(0,2)
ter_edge=zeros(0,2)
non_edge=zeros(0,2)

for i in 1:length(edge_list[:,1])
    if ss_contact_matrix[edge_list[i,1],edge_list[i,2]]!=0
        bar2D+=1
        barTot+=1
        sec_edge=vcat(sec_edge,edge_list[i,:]')
         proxy_4=push!(proxy_4,bar2D)
        proxy_3=push!(proxy_3,bar2D+bar3D)
        proxy_2=push!(proxy_2,bar2D+bar3D+barNeigh)
        proxy_1=push!(proxy_1,barTot)
        
    elseif tertiary_contact_matrix[edge_list[i,1],edge_list[i,2]]!=0
                    bar3D+=1
                    barTot+=1
                    ter_edge=vcat(ter_edge,edge_list[i,:]')
                        proxy_4=push!(proxy_4,bar2D)
                        proxy_3=push!(proxy_3,bar2D+bar3D)
                        proxy_2=push!(proxy_2,bar2D+bar3D+barNeigh)
                        proxy_1=push!(proxy_1,barTot)
              
    elseif edge_list[i,2]-edge_list[i,1]<=4
         barNeigh+=1
        barTot+=1
        
         proxy_4=push!(proxy_4,bar2D)
        proxy_3=push!(proxy_3,bar2D+bar3D)
        proxy_2=push!(proxy_2,bar2D+bar3D+barNeigh)
        proxy_1=push!(proxy_1,barTot)
    else
        barTot+=1
        non_edge=vcat(non_edge,edge_list[i,:]')
         proxy_4=push!(proxy_4,bar2D)
        proxy_3=push!(proxy_3,bar2D+bar3D)
        proxy_2=push!(proxy_2,bar2D+bar3D+barNeigh)
        proxy_1=push!(proxy_1,barTot)
    end
    
end
Plot1=Plots.bar(proxy_1,label="", linewidth=5, linecolor="black", color="black")
Plots.bar!(proxy_1,legend=:topleft,label="OTHER", linewidth=1.1, linecolor="black", color="black")
Plots.bar!(proxy_2,label="NEIGHBOURS", linewidth=1.3, linecolor="lightblue", color="lightblue")
Plots.bar!(proxy_3,label="3D CONTACT",linewidth=1.3, linecolor="ivory", color="ivory")
Plots.bar!(proxy_4,label="SECONDARY STRUCTURE", linewidth=1.3, linecolor="lawngreen", color="lawngreen")

L=len
  
Plot2=Plots.scatter(tertiary_contact_list[1:end,1],tertiary_contact_list[1:end,2],c="white", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4, markerstrokewidth=6,framestyle = :box)
Plots.scatter!(tertiary_contact_list[1:end,1],tertiary_contact_list[1:end,2],c="white", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4.6, markerstrokewidth=0,framestyle = :box)
Plots.scatter!(non_edge[1:end,2],non_edge[1:end,1],color="black", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4.5, markerstrokewidth=6)
Plots.scatter!(ter_edge[1:end,2],ter_edge[1:end,1],color="white", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4, markerstrokewidth=6)
Plots.scatter!(ter_edge[1:end,2],ter_edge[1:end,1],color="white", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4.6, markerstrokewidth=0)
Plots.scatter!(sec_edge[1:end,2],sec_edge[1:end,1],color="white", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4.2, markerstrokewidth=6)
Plots.scatter!(sec_edge[1:end,2],sec_edge[1:end,1],color="green1", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4.6, markerstrokewidth=0)
Plots.scatter!(ss_contact_list[1:end,1],ss_contact_list[1:end,2],c="white", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4.2, markerstrokewidth=6)
Plots.scatter!(ss_contact_list[1:end,1],ss_contact_list[1:end,2],c="green1", aspect_ratio = 1.,xlim = [0, L],
    ylim = [0,  L],legend=false, markersize=4.6, markerstrokewidth=0)
Plots.plot!(LinRange(1,L,L),width=3.2,color="black")
 
return Plots.plot(Plot1,Plot2,layout=(1,2), size=(600,300),legendfont=font(7))
end


