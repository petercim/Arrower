# Arrower
Usage: python arrower_michael.py <GenBank File>\n
  Additional Flags (all in pixels units): 
  -H <ArrowHeight> -E <HeadEdge> -l <HeadLength> -X <marginX> -Y <marginY> -S <scaling> -F <fontSize>\n

  ArrowHeight ... the width of the arrow central part
  HeadEdge    ... additional width of the head
  HeadLength  ... head length
  marginX     ... left-site margins
  marginY     ... top-site margins
  scaling     ... scaling of px per bp (100 means 100bp/px)
  fontSize    ... gene annotation font size\n
  
  After running the script, the SVG script will be printed on the screen. To save
  it directly into the file, one can use command: 
  python arrower_michael.py <file> <flags> > output_file.svg

File can be then further edited in Adobe Illustrator.         
