#would attempt this with pygraphviz but cant fit homebrew on my machine
#as a result the layout could be a little messy
import networkx as nx #use pip install
import matplotlib.pyplot as plt
import pylab


def crudegraph(partlist,linkerlist,linkerbutton):
    
    G = nx.DiGraph()
    linkerlist=['lk0','lk1','lk2','lk3','lk4'] #sample linker parts list
    if linkerbutton==1:
        for n in range ((len(partlist))-1):
            for m in range(len(partlist[n])):
                if n==1:
                    G.add_edges_from([(partlist[n-1][0],partlist[n][m])], weight=1)
                    G.add_edges_from([(partlist[n][m],linkerlist[n])], weight=1)
                elif n>1 and n<len(partlist)-2:
                    G.add_edges_from([(linkerlist[n-1],partlist[n][m])], weight=1)
                    G.add_edges_from([(partlist[n][m],linkerlist[n])], weight=1)
                elif n==len(partlist)-2:
                    G.add_edges_from([(linkerlist[n-1],partlist[n][m])], weight=1)
                    G.add_edges_from([(partlist[n][m],partlist[n+1][0])], weight=1)
    else:
        for n in range (1,(len(partlist))):
            for m in range(len(partlist[n])):
                for o in range(len(partlist[n-1])):
                    G.add_edges_from([(partlist[n-1][o],partlist[n][m])], weight=1)
    
    
    val_map={}            
    for n in range(len(partlist)):
        for m in range(len(partlist[n])):
            val_map.update({partlist[n][m]:n})
    
    values = [val_map.get(node, 0.45) for node in G.nodes()]
    
    pos=nx.spring_layout(G)
    
    labels={}
    
    for node in G.nodes(): 
        labels[node] = node
    
    #ax=plt.axes() 
    plt.figure(figsize=(8,4)) 
    nx.draw(G,pos)        
    nx.draw_networkx_labels(G,pos,labels)
    nx.draw(G,pos, node_color = values, node_size=1500,node_shape='s',edge_cmap=plt.cm.Reds)

    pylab.show()

#crudegraph(partlist,linkerlist,linkerbutton)