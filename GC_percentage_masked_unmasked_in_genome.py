import argparse
import Bio
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter
from matplotlib import font_manager, cm
plt.rcParams["figure.figsize"] = [20.00, 20.00]
import pandas as pd
import numpy as np
import os
import csv
from collections import OrderedDict
from PIL import Image
from glob import glob

def create_output_directory(directory, species, window):
    output_directory = os.path.join(directory, f'{species}_output_{window}')
    os.makedirs(output_directory, exist_ok=True)
    return output_directory

def main(fasta_file, species, window):
    output_directory = create_output_directory(os.path.dirname(fasta_file), species, window)
    print("Output Directory:", output_directory)
    print("Fasta File:", fasta_file)
    print("Species:", species)
    print("Window:", window)
    return output_directory

def how_many(my_seq, string):
    my_sum = 0
    for base in string:
        my_sum += my_seq.count(base)
    return my_sum

# Create an ArgumentParser to handle command line arguments
parser = argparse.ArgumentParser(description="Calculate GC% content of masked and unmasked genome portion")
parser.add_argument("fasta_file", help="Path to the input FASTA file")
parser.add_argument("species", help="Species name")
parser.add_argument("window", type=int, help="Window size for GC% calculation")
args = parser.parse_args()

# Create a list of arguments from the command line arguments
what_to_operate = [(args.fasta_file, args.species, [args.window])]
column_name = 'GC DNA' # replace with your column name
# Set a special value that should be colored differently
special_value = 0.0

bottom = mpl.colormaps['RdYlGn'] #changed from cm.get_cmap('RdYlGn', 255) that is deprectaed from Matplotlib ver 3.7 onwards
special_color = (0.0,0.0,0.999,0.999)
newcolors = np.vstack((special_color,bottom(np.linspace(0, 1, 255))))
new_cmap = ListedColormap(newcolors, name='new_RdYlGn')

slight_move=0.009         

for my_fasta_file in what_to_operate:  #highest level
    for my_window in my_fasta_file[2]:
        print(my_fasta_file[0],my_window)
        sec_values = '' # replace with your release number
        c=[]

        GC_chromosoms=OrderedDict()
        chomosome_no=0

        gc_values=tuple()

        d_sec_values = dict()

        lengths=[]
       

        for rec in SeqIO.parse(f"{my_fasta_file[0]}", "fasta"): 
            GC_chromosoms_windowed={}  # per chromosome
            # count GC values via sliding windows
            c.append(rec.id)
            sec_values += rec.seq
            if len(sec_values) >= 10*my_window:
                gc_values = gc_values + tuple(gc_fraction(sec_values[i:i+my_window]) for i in range(0,len(sec_values),my_window)) # replace with your release number
                sec_values = ''
        
                # chromosomes
                soft_mask= how_many(rec.seq,'acgt')
                ALL_all=how_many(rec.seq,'acgtACGT')
                GC_fraction=gc_fraction(rec.seq)
                chomosome_no += 1
                try:
                    GC_chromosoms[str(chomosome_no) + '+' + rec.id]=((GC_fraction,soft_mask/ALL_all),len(rec.seq))  # we distract N and n   
                except ZeroDivisionError:
                    GC_chromosoms[str(chomosome_no)+'+'+ rec.id]=((GC_fraction,0.0), len(rec.seq)) # distract N and n
                
                # count bases
                for base in 'acgtACGTNn':
                    if base in d_sec_values:
                        d_sec_values[base] += rec.seq.count(base)
                    else:
                        d_sec_values[base] = rec.seq.count(base)

                # max length
                lengths.append(len(rec.seq))
                max_len=max(lengths)
                
                # graphs
                for w in range(0,len(rec.seq),my_window):
              
                    #print(range(i,i+window),len(gc_values_94[i:i+window]))
                    soft_mask= rec.seq[w:w+my_window].count('a') + rec.seq[w:w+my_window].count('t') + rec.seq[w:w+my_window].count('c') + rec.seq[w:w+my_window].count('g') 
                    ALL_all=soft_mask + rec.seq[w:w+my_window].count('A') + rec.seq[w:w+my_window].count('T') + rec.seq[w:w+my_window].count('C') + rec.seq[w:w+my_window].count('G')
                    GC_fraction=gc_fraction(rec.seq[w:w+my_window])
            
                    try:
                        GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,soft_mask/ALL_all)  # distract N and n    
                        #print(range(i,i+window),len(gc_values_94[i:i+window])) 
                    except ZeroDivisionError:
                        GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,0.0)  # distract N and n               
            
                #fig = plt.figure(figsize=(int(50* (len(rec.seq)/max_len)),int(20* (len(rec.seq)/max_len))))
                fig = plt.figure(figsize=(50,20))
                ax1 = fig.add_subplot(111)
                ax1.use_sticky_edges = True
                names = list(GC_chromosoms_windowed.keys())
                v = list(GC_chromosoms_windowed.values())
                v1=[100 * i[0] for i in v]
                #v2=[i[1] for i in v]
                names_part=[n.split('|')[0] for n in names]             
                #values = np.array([special_value if (i[0]+i[1] == 0 and ALL_all== 0) else i[1]+slight_move if (i[1] => 0.0 and i[1] <= 0.00390625 and ALL_all > 0) else i[1] for i in v])  # values for special_indices
                values = np.array([special_value if i[0]+i[1] == 0 else i[1]+slight_move if (i[1] ==0 and i[0] != 0) else i[1] for i in v])  # values for special_indices
       
                if ALL_all== 0:
                    print('blue',rec.id,rec.seq)
                
                sc=ax1.scatter(names_part,v1, s=5, c=values,  cmap=new_cmap, marker="o", label='Blue dots are "NNnn" or spaces' ) # only GC fraction in y axe, only chromosomes
                ax1.grid(True)        
                #plt.xlim(0.0, max_len)
                ax1.xaxis.set_ticks(np.arange(0, max_len//my_window, max_len//12000))  
        
                current_values = plt.gca().get_xticks()
                plt.gca().set_xticklabels(['{:,.0f}'.format(x*1000) for x in current_values])
        
                #ax1.xaxis.set_major_formatter(StrMethodFormatter('{x:,}'))
                ax1.set_facecolor("lightgrey")   
        
                plt.title(f'GC% values of {my_fasta_file[1]} with {my_window} - chromosome {rec.id}', fontsize = 35) # replace the animal name and release no.
                plt.ylabel('GC fraction', fontsize = 30)
                plt.xlabel(f'Chromosome {rec.id} windows', fontsize = 30)
                #plt.xticks(np.arange(0, max_len, max_len//divide_by))
                plt.xticks(fontsize = 20, rotation = 30)
                plt.yticks(fontsize = 30)
                plt.legend(labelcolor='blue')
                
                cbar=plt.colorbar(sc,label="Red 0% soft-masked to green 100% soft-masked", orientation="horizontal")
        
                text = cbar.ax.xaxis.label
                font = font_manager.FontProperties(size=35)
                text.set_font_properties(font)        
        
                tick_font_size = 35
                cbar.ax.tick_params(labelsize=tick_font_size)

                output_directory = main(my_fasta_file[0],my_fasta_file[1],my_window)
                plt.savefig(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_{my_window}_profile_soft_unmask_{rec.id}.png')  # replace the animal name and release no.

                plt.close('all')
                
        columns_dna = [column_name] # replace with your name e.g. release number - default GC_DNA
        gc_df_dna = pd.DataFrame(gc_values, columns = columns_dna) # replace with your release number

        GC_chromosoms = OrderedDict(sorted(GC_chromosoms.items(), key=lambda x: x[1][1],reverse=True))

        #max_len=max(lengths)
    
        gc_df_dna.to_csv(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_per_windows_{my_window}.csv')
    
        #GC_chromosoms # replace with your release number
    
        with open(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_{my_window}.csv','w') as f:  # name csv file as you wish 
            for k,v in GC_chromosoms.items():
                f.write(f'{k.split("+")[0]},{k.split("+")[1]},{v[0][0]},{v[0][1]},{v[1]}\n')
                
        # scattered graphs - color by value, no display
        
        fig = plt.figure(figsize=(15,10))
        ax1 = fig.add_subplot(111)
        ax1.set_facecolor("lightgrey")  
        ax1.use_sticky_edges = True
        #ax1.margins(x=0.0,y=0.0)
        #print(range(i,i+window),len(gc_values_94[i:i+window]))
        
        names = list(i.split('+')[0] for i in GC_chromosoms.keys())
        v = list(GC_chromosoms.values())
        v1=[100 * i[0][0] for i in v]
        v2=[i[0][1] for i in v]

        sc=ax1.scatter(names,v1, s=150, c=v2, cmap='RdYlGn', marker="o", ) # only GC fraction in y axe, only chromosomes
        plt.grid(True)
        plt.title(f'GC% means of size-sorted chromosomes with {my_window} in {my_fasta_file[1]} ', fontsize = 18)  # Replace the animal name
        plt.ylabel('GC fraction', fontsize = 15)
        plt.xlabel('Chromosomes', fontsize = 15)
        plt.xticks(rotation = 30, fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.colorbar(sc,label="Red 0% soft-masked to green 100% soft-masked", orientation="horizontal")
        plt.savefig(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_soft_unmask_per_chromosomes_{my_window}_{my_fasta_file[1]}.png')  # Replace the animal name
        #plt.show()
        plt.close('all')
        #print(f'Done {my_fasta_file[1]}')
        
        #3EVAN
        # d sec values
        with open(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_{my_window}_something.csv','w') as f:  # name csv file as you wish 
            for k,v in d_sec_values.items():
               f.write(f'{k},{v}\n')
                
        # results % GC
        gc_all = (d_sec_values['c']+d_sec_values['g']+d_sec_values['C']+d_sec_values['G'])/(d_sec_values['c']+d_sec_values['g']+d_sec_values['C']+d_sec_values['G']+ d_sec_values['a']+d_sec_values['t']+d_sec_values['A']+d_sec_values['T'])
       
        with open(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_all_{my_window}_GCproc.csv','w') as f:  # name csv file as you wish 
            f.write(f'%GC of {my_fasta_file[1]} is {gc_all}\n')
            
        #  Soft vs Unmasked
        # GC to ATGC, gc to atgc
        GC_to_ATGC={}
        atgc_lower= d_sec_values['a'] + d_sec_values['t'] + d_sec_values['c']  + d_sec_values['g']  # 'atgc'
        gc_lower=d_sec_values['c']  + d_sec_values['g'] # cg
        atgc_upper=d_sec_values['A'] + d_sec_values['T'] + d_sec_values['C']  + d_sec_values['G'] # ATCG
        gc_upper=d_sec_values['C']  + d_sec_values['G'] # CG
        GC_to_ATGC['lower']=gc_lower/atgc_lower  # distract N and n
        GC_to_ATGC['upper']=gc_upper/atgc_upper  # distract N and n
        GC_to_ATGC['repeat%']=atgc_lower/(atgc_lower + atgc_upper) #calculated the masked/repetitive fraction   
        GC_to_ATGC.items()
        
       	with open(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_lower_upper_repeat_{my_window}.csv','w') as f:  # name csv file as you wish 
            for k,v in GC_to_ATGC.items():
               f.write(f'{k},{v}\n')
                
        # GC soft and (small) and unmasked (capital) repetitions - by chromosomes
        # GC DNA release ..... - by chromosomes

# adjusted chromosome density
for my_fasta_file in what_to_operate:  #highest level
    for my_window in my_fasta_file[2]:
        print(my_fasta_file[0],my_window)
        sec_values = '' # replace with your release number
        c=[]

        GC_chromosoms=OrderedDict()
        chomosome_no=0

        gc_values=tuple()

        d_sec_values = dict()

        lengths=[]
		
        for rec in SeqIO.parse(f"{my_fasta_file[0]}", "fasta"): 
            GC_chromosoms_windowed={}  # per chromosome
            #if eval(my_fasta_file[3]):  # chromosomes only
            # count GC values via sliding windows
            c.append(rec.id)
            sec_values += rec.seq
            if len(sec_values) >= 10*my_window:
                gc_values = gc_values + tuple(gc_fraction(sec_values[i:i+my_window]) for i in range(0,len(sec_values),my_window)) # replace with your release number
                sec_values = ''
        
            # chromosomes
            soft_mask= how_many(rec.seq,'acgt')
            ALL_all=how_many(rec.seq,'acgtACGT')
            GC_fraction=gc_fraction(rec.seq)
            chomosome_no += 1
                
            try:
                GC_chromosoms[str(chomosome_no) + '+' + rec.id]=((GC_fraction,soft_mask/ALL_all),len(rec.seq))  # we distract N and n   
            except ZeroDivisionError:
                GC_chromosoms[str(chomosome_no)+'+'+ rec.id]=((GC_fraction,0.0), len(rec.seq)) # distract N and n
                  
        
                # count bases
            for base in 'acgtACGTNn':
                if base in d_sec_values:
                    d_sec_values[base] += rec.seq.count(base)
                else:
                    d_sec_values[base] = rec.seq.count(base)

            # max length
            lengths.append(len(rec.seq))
            max_len=max(lengths)
                
                #print(rec.id)  
            for w in range(0,len(rec.seq),my_window):
              
                    #print(range(i,i+window),len(gc_values_94[i:i+window]))
                soft_mask=sum(1 for i in rec.seq[w:w+my_window] if i in 'acgt')
                ALL_all=sum(1 for i in rec.seq[w:w+my_window] if i in 'acgtACGT')
                GC_fraction=gc_fraction(rec.seq[w:w+my_window])
            
                try:
                    GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,soft_mask/ALL_all)  # distract N and n    
                        #print(range(i,i+window),len(gc_values_94[i:i+window])) 
                except ZeroDivisionError:
                    GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,0.0)  # distract N and n
            
                #fig = plt.figure(figsize=(int(50* (len(rec.seq)/max_len)),int(20* (len(rec.seq)/max_len))))
            fig = plt.figure(figsize=(25,10))
            ax1 = fig.add_subplot(111)
        
            names = list(GC_chromosoms_windowed.keys())
            v = list(GC_chromosoms_windowed.values())
            v1=[100 * i[0] for i in v]
                #v2=[i[1] for i in v]
            names_part=[n.split('|')[0] for n in names]
            values = np.array([special_value if i[0]+i[1] == 0 else i[1]+slight_move if (i[1] ==0 and i[0] != 0) else i[1] for i in v])  # values for special_indices
       
            sc=ax1.scatter(names_part,v1, s=5, c=values, cmap=new_cmap, marker="o", label='Blue dots are "NNnn" or spaces' ) # only GC fraction in y axe, only chromosomes
                #plt.xlim(0.0, max_len)
            ax1.xaxis.set_ticks(np.arange(0, max_len//my_window, max_len//12000))    
            plt.title(f'GC% values of {my_fasta_file[1]} - {my_window} chromosome {rec.id}', fontsize = 30) # replace the animal name and release no.  
            plt.xticks(fontsize = 20, rotation = 30)
            plt.yticks(fontsize = 20)
            plt.legend(labelcolor='blue')
            plt.savefig(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_{my_window}_nobar_soft_unmask_{rec.id}.{len(rec.seq)}.png')  # replace the animal name and release no.
			
                #plt.show()
            plt.close('all')

output_directory = main(my_fasta_file[0],my_fasta_file[1],my_window)

def get_concat(all_images,im_width,im_height):
    
    images = [Image.open(all_images[i]) for i in range(len(all_images))]
    if len(all_images)%2 == 0:
        how_high=len(all_images)
    else:
        how_high=len(all_images)+1
        
    dst = Image.new('RGB', (2*im_width, int((how_high * im_height)/2)))
    dst.paste(Image.open(f'{output_directory}/{os.path.basename(my_fasta_file[1])}headtitle_{my_fasta_file[1]}_{my_window}.png'), (0, 0))
        
    for i in range(0,len(images)+1,2):
        try:
            dst.paste(images[i], (0, int((i * im_height)/2)))
        except IndexError:
            pass
        try:            
            dst.paste(images[i+1], (im_width,  int((i * im_height)/2)))
        except IndexError:
            pass
    return dst

for my_fasta_file in what_to_operate:  #highest level
    for my_window in my_fasta_file[2]:
        #print(directories_only(my_fasta_file[0]),my_window)
        # cycle for all nobar images for species and window
        all_images = sorted([f for f in glob(f"{output_directory}/{os.path.basename(my_fasta_file[1])}*{my_window}_nobar*.png")],key=lambda x: os.path.getsize(x), reverse = True)
        
        images = [Image.open(f"{x}") for x in all_images]
        widths, heights = zip(*(i.size for i in images))
        max_width = max(widths)
        total_height = sum(heights)
        new_im = Image.new('RGB', (max_width, total_height))

        y_offset = 0
        for im in images:
            new_im.paste(im, (0,y_offset,))
            y_offset += im.size[1]

        new_im.save(f'{output_directory}/{os.path.basename(my_fasta_file[1])}{my_fasta_file[1]}_{my_window}_all_graphs.png')
        
        # Read small graphs one by one and merge them all to two columns¶
        #fig = plt.figure(figsize=(0.01,0.01))
        #plt.title(f'{my_fasta_file[1]} - sliding window {my_window} ', fontsize = 18)  # Replace the animal name

        #plt.savefig(f'{output_directory}/{os.path.basename(my_fasta_file[1])}_headtitle_{my_window}.png')  # Replace the animal name
        #plt.show()
        
        #plt.close('all')
        
        im = Image.open(all_images[0])
        im_width=im.size[0]
        im_height=im.size[1]
        im.close()
        
        #output_directory = main(my_fasta_file[0],my_fasta_file[1],my_window)
        #output_directory = main(args.fasta_file, args.species, [args.window])
		#os.path.dirname(output_directory)

        def save_concatenated_images(output_directory, filename, all_images):
        # Set the number of columns and rows for the grid
            columns = 2
            rows = (len(all_images) + columns - 1) // columns

            # Create a new figure with subplots
            fig, axes = plt.subplots(rows, columns, figsize=(100,100))

            # Flatten the axes array to handle cases with a single row or column
            axes = axes.flatten()

            # Loop through the images and plot them in the subplots
            for i, image_path in enumerate(all_images):
                if i < len(axes):
                    img = plt.imread(image_path)
                    axes[i].imshow(img)
                    axes[i].axis('off')

            # Adjust layout and save the figure
            plt.subplots_adjust(wspace=0, hspace=0)
            
            plt.savefig(os.path.join(output_directory, filename), bbox_inches='tight', pad_inches=0)
            
            plt.close()
			
            #plt.tight_layout()
            #plt.savefig(os.path.join(output_directory, filename))

        # Example usage
        output_directory = main(my_fasta_file[0],my_fasta_file[1],my_window)
        filename = f'{my_fasta_file[1]}_all_small_graphs_{my_window}.png'
        save_concatenated_images(output_directory, filename, all_images)

        #get_concat(all_images, im_width, im_height).save(os.path.join(output_directory, f'{my_fasta_file[1]}_all_small_graphs_{my_window}.png'))

#get_concat(all_images,im_width,im_height).save(f'{os.path.dirname(output_directory)(my_fasta_file[1])}_all_small_graphs_{my_window}.png')
        #print(f'{my_window}_{my_fasta_file[1]}')
	
if __name__ == "__main__":
    #parser.add_argument("arguments", type=str, help="Provide 'fasta_file_path,species,window' as a single string argument")
    output_directory = create_output_directory(os.path.dirname(args.fasta_file), args.species, args.window)
    main(args.fasta_file, args.species, [args.window])
    
