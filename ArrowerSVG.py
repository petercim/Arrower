#! /usr/bin/python

######################################################################
#                                                                    #
#           PLOT ARROWS FOR GENE CLUSTER GIVEN A GenBank FILE        #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                                                                    #
######################################################################

from Bio import SeqIO
import sys

# --- Draw arrow for gene
def arrow(X,Y,L,H,strand,h,l,color):
    '''
    SVG code for arrow:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color
        - strand
    the edges are ABCDEFG starting from (X,Y)     
    '''
    
    if strand == '+':
        
        A = [X,Y]
        B = [X+L-l,Y]
        C = [X+L-l,Y-h]
        D = [X+L,Y+H/2]
        E = [X+L-l,Y+H+h]
        F = [X+L-l,Y+H]
        G = [X,Y+H]

        if L < l:
            # squize arrow if length shorter than head length
            B = [X,Y]
            C = [X,Y-h]
            D = [X+L,Y+H/2]
            E = [X,Y+H+h]
            F = [X,Y+H]

        line = """ <polygon id="bla" points="%i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i"
                   style="fill:#cccccc;fill-opacity:1.0;
                   stroke:#000000;stroke-width:2">
                   </polygon>""" % (A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],E[0],E[1],F[0],F[1],G[0],G[1])
        
        return line

    elif strand == '-':
        
        A = [X+L,Y]
        B = [X+l,Y]
        C = [X+l,Y-h]
        D = [X,Y+H/2]
        E = [X+l,Y+H+h]
        F = [X+l,Y+H]
        G = [X+L,Y+H]

        if L < l:
            # squize arrow if length shorter than head length
            B = [X+L,Y]
            C = [X+L,Y-h]
            D = [X,Y+H/2]
            E = [X+L,Y+H+h]
            F = [X+L,Y+H]

        line = """ <polygon id="bla" points="%i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i"
                   style="fill:#cccccc;fill-opacity:1.0;
                   stroke:#000000;stroke-width:2">
                   </polygon>""" % (A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],E[0],E[1],F[0],F[1],G[0],G[1])
        
        return line
    
    else: return 0

def line(X,Y,L):
    '''
    Draw a line below genes
    '''
    
    line = """<line x1="%i" y1="%i" x2="%i" y2="%i"
              style="stroke:rgb(99,99,99);stroke-width:2"/>""" % (X,Y,X+L,Y)
    return line


def text(text,X,Y,F):
    '''
    Write gene name
    '''
    line = """<text x="%i"  y="%i"
          style="font-family: Arial;
                 font-size  : %i;
                 font-style : italic;
                "
          >%s</text>""" % (X,Y,F,text)
    return line
    

def RGBToHTMLColor(rgb_tuple):
    '''
    Convert RGB color to HTML color format
    '''
    hexcolor = '#%02x%02x%02x' % rgb_tuple
    
    return hexcolor.strip('-')


def SVG(GenBankFile,ArrowHeight=20,HeadEdge=8,HeadLength=10,marginX=100,marginY=30,scaling=100.0,font=14):
    '''
    Create the main SVG document:
        - read in GenBank documnet
        - find genes, start and stop positions, and strands
        - write the SVG files
    '''
    
    # --- create SVG header
    ALL_TEXT = ""

    header = """<?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

    <svg width="870" height="300" xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'
    onload='Init(evt)'>
    """ #% Width

    ALL_TEXT += header

    # --- read in GenBank file
    file = open(GenBankFile,'r')
    for seq_record in SeqIO.parse(file, "genbank"):
        
        # draw a line that corespond to cluster size
        ClusterSize = len(seq_record.seq)
        ALL_TEXT += line(marginX,marginY+ArrowHeight/2,ClusterSize / scaling)
        
        previousStop = 0 # keep track of previous stop to adjust text location
        previousLevel = 1
        for feature in seq_record.features:
            if feature.type == 'CDS':
                
                GeneName = feature.qualifiers['gene'][0]

                strand = feature.strand
                if strand == -1:
                    strand = '-'
                elif strand == 1:
                    strand = '+'

                start = str(feature.location.start)
                if '>' in start or '<' in start:
                    start = start.replace('>','').replace('<','')
                start = int(start) / scaling
                
                stop = str(feature.location.end) 
                if '>' in stop or '<' in stop:
                    stop = stop.replace('>','').replace('<','')
                stop = int(stop) / scaling
                
                # write arrow to SVG file
                ALL_TEXT += arrow(start+marginX,marginY,stop-start,ArrowHeight,strand,HeadEdge,HeadLength,GeneName)

                # annotate genes
                if previousLevel  == 4: previousLevel = 0
                if previousStop == 0:
                    ALL_TEXT += text(GeneName,marginX+start+(stop-start)/2-len(GeneName)*F/4.,marginY+50,F)
                else:
                    if ((start+(stop-start)/2-len(GeneName)*F/4.) - previousStop) < 10:
                        ALL_TEXT += text(GeneName,marginX+start+(stop-start)/2-len(GeneName)*F/4.,marginY+50+(F+6)*previousLevel,F)
                        previousLevel += 1
                    else:
                        previousLevel = 0
                        ALL_TEXT += text(GeneName,marginX+start+(stop-start)/2-len(GeneName)*F/4.,marginY+50+(F+6)*previousLevel,F)
                previousStop = stop


    ALL_TEXT += '</svg>'
   
    print ALL_TEXT

ARGS = sys.argv
try:
    if len(ARGS) < 2:
        print '''
          Correct usage: python arrower_michael.py <GenBank File>\n
          \tAdditional Flags (all in pixels units): 
          \t  -H <ArrowHeight> -E <HeadEdge> -l <HeadLength> -X <marginX> -Y <marginY> -S <scaling> -F <fontSize>\n
          \t ArrowHeight ... the width of the arrow central part
          \t HeadEdge    ... additional width of the head
          \t HeadLength  ... head length
          \t marginX     ... left-site margins
          \t marginY     ... top-site margins
          \t scaling     ... scaling of px per bp (100 means 100bp/px)
          \t fontSize    ... gene annotation font size\n
          After running the script, the SVG script will be printed on the screen. To save
          it directly into the file, one can use command: 
          \tpython arrower_michael.py <file> <flags> > output_file.svg
          File can be then further edited in Adobe Illustrator.         
          '''

    else:
        file = ARGS[1]  # GenBank file
        H = 20          # arrow width
        E = 8           # arrow head edge size
        l = 10          # arrow head edge length
        X = 100         # left-site margins
        Y = 30          # top-site margins
        S = 100.0       # scaling factor
        F = 14

        for x in xrange(len(ARGS)):
            if '-' in ARGS[x]:
                if 'H' in ARGS[x]: H = int(ARGS[x+1])
                elif 'E' in ARGS[x]: E = int(ARGS[x+1])
                elif 'l' in ARGS[x]: l = int(ARGS[x+1])
                elif 'X' in ARGS[x]: X = int(ARGS[x+1])
                elif 'Y' in ARGS[x]: Y = int(ARGS[x+1])
                elif 'S' in ARGS[x]: S = float(ARGS[x+1])
                elif 'F' in ARGS[x]: F = int(ARGS[x+1])
                else:
                    print 'Check your flags!'
                    sys.exit(1)
        SVG(file,ArrowHeight=H,HeadEdge=E,HeadLength=l,marginX=X,marginY=Y,scaling=S,font=F)

except:
    print ''' Incorrect usage. Please, read manual again:
          Correct usage: python arrower_michael.py <GenBank File>\n
          \tAdditional Flags (all in pixels units): 
          \t  -H <ArrowHeight> -E <HeadEdge> -l <HeadLength> -X <marginX> -Y <marginY> -S <scaling> -F <fontSize>\n
          \t ArrowHeight ... the width of the arrow central part
          \t HeadEdge    ... additional width of the head
          \t HeadLength  ... head length
          \t marginX     ... left-site margins
          \t marginY     ... top-site margins
          \t scaling     ... scaling of px per bp (100 means 100bp/px)
          \t fontSize    ... gene annotation font size\n
          After running the script, the SVG script will be printed on the screen. To save
          it directly into the file, one can use command: 
          \tpython arrower_michael.py <file> <flags> > output_file.svg
          File can be then further edited in Adobe Illustrator.         
          ''' 
    sys.exit(1)
