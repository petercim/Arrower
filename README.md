# Arrower
Usage: python arrower_michael.py <GenBank File>

  Additional Flags (all in pixels units): 

  - -H: ArrowHeight ... the width of the arrow central part
  - -E: HeadEdge    ... additional width of the head
  - -l: HeadLength  ... head length
  - -X: marginX     ... left-site margins
  - -Y: marginY     ... top-site margins
  - -S: scaling     ... scaling of px per bp (100 means 100bp/px)
  - -F: fontSize    ... gene annotation font size\n
  
After running the script, the SVG script will be printed on the screen. To save it directly into the file, one can use command: 
  
 - python arrower_michael.py <file> <flags> > output_file.svg

File can be then further edited in Adobe Illustrator.         
