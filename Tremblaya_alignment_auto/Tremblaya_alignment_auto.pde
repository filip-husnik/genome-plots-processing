import processing.svg.*;



//Input is a CSV file called "data_files.csv" placed in the data folder where:
//                        column 1 = genome sizes in base pairs, sorted on the host phylogeny (genome)
//                        column 2 = path to files with gene coordinates and functions (genes)
//                        column 3 = path to files with alignment connections (connections)


// Other global variables
Table gene_table;
float figure_width = 1800; // 100 pixels smaller than the canvas width to have some free space
float figure_height = 1000; // same as canvas height
int height_unit = 5;
int gene_height = 3 * height_unit;
int genome_height = 10 * height_unit;



void setup() {
  noLoop();
  size(1900, 1000); // see above
  background(255);
  strokeWeight(0.003);
  beginRecord(SVG, "genome_alignment.svg");
  gene_table = loadTable("./data/data_files.csv", "header");
}



void drawGenes(String FileName, float genome_size, int y) {

  int[] genome = gene_table.getIntColumn(0);
  float max_genome = max(genome);
  float ratio = max_genome/figure_width; //(bps per pixel)
  
  float offset = figure_height-(figure_height/(max_genome/genome_size));
  Table table = loadTable(FileName);
  int numEntries = table.getRowCount();
  for (int a = 0; a < numEntries; a += 1) {
    float start_coord = table.getFloat(a, 0);
    float stop_coord = table.getFloat(a, 1);
    String gene_type = table.getString(a, 2);
    String strand = table.getString(a, 3);
    float x1 = start_coord/ratio;
    float x2 = stop_coord/ratio;
    float x3 = x2-x1; //gene length
    
    if (gene_type.equals("CDS")) {
      fill(0, 0, 255);
    } else if (gene_type.equals("pseudogene")) {
      fill(64, 64, 64);
    } else if (gene_type.equals("mature_transcript") || gene_type.equals("misc_RNA")) {
      fill(255, 255, 0);
    } else if (gene_type.equals("tRNA")) {
      fill(255, 0, 0);
    } else if (gene_type.equals("rRNA")) {
      fill(173, 255, 47);
    }
    if (strand.equals("-")) {
      rect(x1+(genome_height+offset), y-10+gene_height, x3, gene_height);
    } else if (strand.equals("+")) {  
      rect(x1+(genome_height+offset), y-5-gene_height, x3, gene_height);
    }
     //genome line
    line(offset+genome_height, y, (genome_size/ratio)+offset+genome_height, y);
    noFill();
    stroke(153, 30);
    //genome box
    rect(offset+genome_height, y-(genome_height/2), genome_size/ratio, genome_height);
  }
}

void drawConnections(String FileName, int length1, int length2, int y1, int y2) {
  
  int[] genome = gene_table.getIntColumn(0);
  float max_genome = max(genome);
  float ratio = max_genome/figure_width; //(bps per pixel)

  Table table = loadTable(FileName);
  print(FileName + "\n");
  int numEntries = table.getRowCount();
  for (int a = 0; a < numEntries; a += 1) {
    float start_coord1 = table.getFloat(a, 0);
    float stop_coord1 = table.getFloat(a, 1);
    float start_coord2 = table.getFloat(a, 3);
    float stop_coord2 = table.getFloat(a, 4);
    String strand = table.getString(a, 2);
    float x1_symb1 = start_coord1/ratio;
    float x2_symb1 = stop_coord1/ratio;
    float x1_symb2 = start_coord2/ratio;
    float x2_symb2 = stop_coord2/ratio;
    float offset1 = figure_height-(figure_height/(max_genome/length1));
    float offset2 = figure_height-(figure_height/(max_genome/length2));


    stroke(255, 0, 0, 20);
    fill(255, 0, 0, 20);
    
    if (strand.equals("-")) {
      quad(x2_symb1+genome_height+offset1, y1, x1_symb1+genome_height+offset1, y1, x1_symb2+genome_height+offset2, y2, x2_symb2+genome_height+offset2, y2);
    } else if (strand.equals("+")) {  
      quad(x2_symb1+genome_height+offset1, y1, x1_symb1+genome_height+offset1, y1, x1_symb2+genome_height+offset2, y2, x2_symb2+genome_height+offset2, y2);
    }
    
  }
}



void draw() {
  
  //To draw connections (if there are any)
  String[] connections = gene_table.getStringColumn(2);
  for (int i = 0; i < connections.length; i++) {
    if ((connections[i] != null) && (connections[i] != "")) {
      drawConnections(gene_table.getString(i,2), gene_table.getInt(i,0), gene_table.getInt(i+1,0), (i+1)*2*genome_height, (i+2)*2*genome_height);
    
    }
  }
  
  //To draw just genes
  int num_genes = gene_table.getIntColumn(0).length;
  if (num_genes == gene_table.getStringColumn(1).length) {
    for (int i = 0; i < num_genes; i++) {
      drawGenes(gene_table.getString(i,1), gene_table.getInt(i,0), (i+1)*2*genome_height);
    }
  }
  

  save("genome_alignment.tif");
  endRecord();
}
