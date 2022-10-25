import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

#dictionary that stores index of key of annotation dataframe
#given a combo of gene+go+evidence -> index of row will be given
annotation_dict = {}

neighbors = {}

#convert a list into a string
#ex. [mitochondrion,localization] = mitochondrion localization
def list_as_string(list):
    string = ""
    for elem in list:
        string += f"{elem} "
    return string[:-1]

#parses the go file and returns a graph wit nodes containing information
def go_parser():
    global neighbors
    g = nx.DiGraph()
    go_file = open("go.obo",'r')
    go_lines = go_file.readlines()
    start = False
    node_id = ""
    node_parts = {}
    node_is = {}
    for line in go_lines:
        if start:
            data = line.split(" ")
            if data[0] == 'id:':
                if data[1] not in g.nodes():
                    node_id = data[1][:-1]
                else:
                    start = False
            elif data[0] == 'is_a:':
                node_is[data[1]] = list_as_string(data[3:])[:-1]
            elif data[0] == 'relationship:' and data[1] == "part_of":
                node_parts[data[2]] = list_as_string(data[4:])[:-1]
        if line == "[Term]\n":
            if start:
                g.add_node(node_id, part_of = node_parts, is_a = node_is)
                node_id = ""
                node_parts = {}
                node_is = {}
            start = True
        elif line[0] == '[':
            start = False
    
    part_of = nx.get_node_attributes(g, "part_of")
    is_a = nx.get_node_attributes(g, "is_a")
    for u in g.nodes():
        if u not in neighbors:
            neighbors[u] = []
        for v in part_of[u]:
            neighbors[u].append(v)
            g.add_edge(u,v)
        for v in is_a[u]:
            neighbors[u].append(v)
            g.add_edge(u,v)
    return g

def annotation_parser():
    global annotation_dict

    go_file = open("goa_human.gaf",'r')
    go_lines = go_file.readlines()
    result = {}

    count = 0
    for line in go_lines:
        if line[0] != '!':
            #1,3,4,6,8 <- [1,n]
            #1,4,6 <- Gene ID, GO ID, Evidence Code
            data = line.split("\t")
            if data[3][0:4] != 'NOT|':
                if f'{data[1]}{data[4]}{data[6]}' not in annotation_dict:
                    if data[4] not in result:
                        result[data[4]] = {}
                    annotation_dict[f'{data[1]}{data[4]}{data[6]}'] = count
                    result[data[4]][data[1]] = [] #gene is key to get other information
                    result[data[4]][data[1]] .append(data[3]) #qualifer = 0
                    result[data[4]][data[1]] .append(data[6]) #evidence = 1
                    result[data[4]][data[1]] .append(data[8]) #aspect = 2
                    count += 1

    #dictionary['go_id']['gene_id'] = [qualifer,evidence,aspect]
    return result

def transfer_up_annotations(go_dag, human_annotations):
    '''
    ancestor_graphs = {}
    processed = {}
    queue = []
    part_of = nx.get_node_attributes(go_dag, "part_of")
    is_a = nx.get_node_attributes(go_dag, "is_a")

    for t in go_dag.nodes():
        processed[t] = 0
        ancestor_graphs[t] = list(nx.ancestors(go_dag, t))
        if len(part_of[t]) == 0 and len(is_a[t]):
            processed[t] = 1
            queue.append(t)
    '''
    #in a sense, since all gene interactions go up the tree, each node must have the same
    #interactions as all there descendants (opposite of ancestors), so we are transfering up
    #gene interactions by making sure that each node has gene interactions of all of it's decendents
    for node in go_dag.nodes():
        if node not in human_annotations:
            human_annotations[node] = {}
        for decendent in list(nx.descendants(go_dag, node)):
            if decendent in human_annotations:
                for gene in human_annotations[decendent]:
                    human_annotations[node][gene] = human_annotations[decendent][gene]
    
    return human_annotations

def get_most_specific(aspect,go_dag, annotations, threshold):
    if aspect == "GO:0008150":
        identifier = 'P'
    elif aspect == "GO:0003674":
        identifier = 'F'
    else:
        identifier = 'C'

    s_subset = []
    minimum = int(len(annotations[aspect])*threshold)
    for t in go_dag:
        count = 0
        if t in annotations:
            for gene in annotations[t]:
                if gene[2] == identifier:
                    count += 1
            if count >= minimum:
                s_subset.append(t)
    t_subset = []
    ancestors= {}
    for t in s_subset:
        ancestors[t] = list(nx.ancestors(go_dag, t))
    
    for t in s_subset:
        check = True
        for d in s_subset:
            if d != t and t in ancestors[d]:
                check = False
                break
        if check:
            t_subset.append(t)
    
    return [s_subset,t_subset] 

def evidence_grouping(annotations):
    exp_codes = ['EXP','IDA','IPI','IMP','IGI','IEP']
    htp_codes = ['HTP', 'HDA', 'HMP', 'HGI', 'HEP']
    phylo_codes = ['IBA', 'IBD', 'IKR', 'IRD']
    comp_codes = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA']
    auth_codes = ['TAS', 'NAS','IC']
    for go_id in annotations:
        for gene in annotations[go_id]:
            evidence_code = annotations[go_id][gene][1]
            if evidence_code in exp_codes:
                annotations[go_id][gene][1] = 'EXP'
            elif evidence_code in htp_codes:
                annotations[go_id][gene][1] = 'HTP'
            elif evidence_code in phylo_codes:
                annotations[go_id][gene][1] = 'PHYLO'
            elif evidence_code in comp_codes:
                annotations[go_id][gene][1] = 'COMP'
            elif evidence_code in auth_codes:
                annotations[go_id][gene][1] = 'AUTH'
            #if the evidence codes are IEA or ND, they remain so (only one evidence code for them anyways)
    return annotations

def gene_evidence_group_counting(outfile, annotations, t_subset):
    outfile.write(f'exp_count,htp_count,phylo_count,comp_count,auth_count,iea_count,nd_count\n')
    exp_count = 0 
    htp_count = 0 
    phylo_count = 0 
    comp_count = 0 
    auth_count = 0 
    iea_count = 0
    nd_count = 0
    assesed = []
    for go_id in t_subset:
        for gene in annotations[go_id]:
            if gene not in assesed:
                if annotations[go_id][gene][1] == 'EXP':
                    exp_count += 1
                elif annotations[go_id][gene][1] == 'HTP':
                    htp_count += 1
                elif annotations[go_id][gene][1] == 'PHYLO':
                    phylo_count += 1
                elif annotations[go_id][gene][1] == 'COMP':
                    comp_count += 1
                elif annotations[go_id][gene][1] == 'AUTH':
                    auth_count += 1
                elif annotations[go_id][gene][1] == 'IEA':
                    iea_count += 1
                elif annotations[go_id][gene][1] == 'ND':
                    nd_count += 1
                assesed.append(gene)
    outfile.write(f'{exp_count},{htp_count},{phylo_count},{comp_count},{auth_count},{iea_count},{nd_count}')

def plot_maker():
    one_percent = ['output_bp_1','output_mf_1','output_cc_1']
    ten_percent = ['output_bp_10','output_mf_10','output_cc_10']
    x = ['Biological Process', 'Molecular Function', 'Cellular Component']

    for threshold in [one_percent,ten_percent]:
        exp= [] 
        htp= [] 
        phylo= [] 
        comp= [] 
        auth= [] 
        iea= []
        nd= []
        for i in threshold:
            outfile = open(f'outputs/{i}', 'r')
            lines = outfile.readlines()
            count = 0
            for line in lines:
                if count == 1:
                    data = line.split(',')
                    exp.append(int(data[0])) 
                    htp.append(int(data[1])) 
                    phylo.append(int(data[2])) 
                    comp.append(int(data[3])) 
                    auth.append(int(data[4])) 
                    iea.append(int(data[5]))
                    nd.append(int(data[6]))
                count += 1
        exp= np.array(exp) 
        htp= np.array(htp) 
        phylo= np.array(phylo) 
        comp= np.array(comp) 
        auth= np.array(auth) 
        iea= np.array(iea)
        nd= np.array(nd)

        plt.bar(x, exp, color='b')
        plt.bar(x, htp, bottom=exp, color='g')
        plt.bar(x, phylo, bottom=exp+htp, color='r')
        plt.bar(x, comp, bottom=exp+htp+phylo, color='c')
        plt.bar(x, auth, bottom=exp+htp+phylo+comp, color='m')
        plt.bar(x, iea, bottom=exp+htp+phylo+comp+auth, color='y')
        plt.bar(x, nd, bottom=exp+htp+phylo+comp+auth+iea, color='k')

        plt.xlabel("Aspect")
        plt.ylabel("Frequency")
        plt.legend(["EXP", "HTP", "PHYLO", "COMP", "AUTH", "IEA", "ND"])
        percent = 1
        if threshold == ten_percent:
            percent = 10
        plt.title(f"Number of Annotations by Evidence ({percent}% Threshold)")
        plt.show()

def main():
    go_dag = go_parser()
    human_annotations = annotation_parser()
    human_annotations_transferred = transfer_up_annotations(go_dag,human_annotations)
    aspect_roots = ["GO:0008150", "GO:0003674", "GO:0005575"] #bp, mf, cc
    most_specific = []
    for threshold in [0.01,0.1]:
        for aspect in aspect_roots:
            most_specific.append(get_most_specific(aspect,go_dag,human_annotations_transferred,threshold))
    human_annotations_transferred = evidence_grouping(human_annotations_transferred)
    names = ['bp_1','mf_1','cc_1','bp_10','mf_10','cc_10']
    i = 0
    for t in most_specific:
        outfile = open(f'output_{names[i]}','w')
        gene_evidence_group_counting(outfile,human_annotations_transferred,t[1])
        i += 1
    plot_maker()

main()