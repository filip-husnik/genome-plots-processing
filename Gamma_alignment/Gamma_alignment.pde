import processing.svg.*;

// Genome sizes and file names in the following three arrays have to be modified based on your data
int[] genome = {
  834734, //DEMHIR
  938041, //GEFVIR
  538294, //MEPCIT
  352837, //MEPMAR
  628221, //HETPER
}; //genome sizes in base pairs, sorted based on the host phylogeny (the same way as files below)

String[] genes = {
  "./data/DEMHIR_gene_info_processing.csv", 
  "./data/GEFVIR_gene_info_processing.csv", 
  "./data/MEPCIT_gene_info_processing.csv", 
  "./data/MEPMAR_gene_info_processing.csv", 
  "./data/HETPER_gene_info_processing.csv", 
}; //gene coordinates and functions

String[] connections = {
  "./data/DEMHIR_vs_GEFVIR.csv", 
  "./data/GEFVIR_vs_MEPCIT.csv", 
  "./data/MEPCIT_vs_MEPMAR.csv", 
  "./data/MEPMAR_vs_HETPER.csv", 
}; //alignment connections

// Other global variables
float max_genome = genome[1];
int genome_nr = genome.length;
float figure_width = 1800; // 100 pixels smaller than the canvas width to have some free space
float figure_height = 1000; // same as canvas height
float ratio = max_genome/figure_width; //(bps per pixel)
int height_unit = 5;
int gene_height = 3 * height_unit;
int genome_height = 10 * height_unit;

void setup() {
  noLoop();
  size(1900, 1000); // see above
  background(255);
  strokeWeight(0.003);
  beginRecord(SVG, "genome_alignment.svg");
}

void drawGenes(String FileName, float genome_size, int y) {
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
      fill(237, 106, 107);
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
  Table table = loadTable(FileName);
  int numEntries = table.getRowCount();
  for (int a = 0; a < numEntries; a += 1) {
    float start_coord1 = table.getFloat(a, 1);
    float stop_coord1 = table.getFloat(a, 2);
    float start_coord2 = table.getFloat(a, 5);
    float stop_coord2 = table.getFloat(a, 6);
    String gene_type1 = table.getString(a, 3);
    String gene_type2 = table.getString(a, 7);
    String strand1 = table.getString(a, 4);
    String strand2 = table.getString(a, 8);
    float x1_symb1 = start_coord1/ratio;
    float x2_symb1 = stop_coord1/ratio;
    float x1_symb2 = start_coord2/ratio;
    float x2_symb2 = stop_coord2/ratio;
    float offset1 = figure_height-(figure_height/(max_genome/length1));
    float offset2 = figure_height-(figure_height/(max_genome/length2));

    if (gene_type1.equals("CDS")) {
      stroke(237, 106, 107, 40);
      fill(237, 106, 107, 40);
    } else if (gene_type1.equals("pseudogene") || gene_type2.equals("pseudogene")) {
      stroke(64, 64, 64, 20);
      fill(64, 64, 64, 20);
    } else if (gene_type1.equals("mature_transcript") || gene_type1.equals("misc_RNA")) {
      stroke(255, 255, 0, 120);
      fill(255, 255, 0, 120);
    } else if (gene_type1.equals("tRNA") || gene_type2.equals("tRNA")) {
      noFill();
    } else if (gene_type1.equals("rRNA")) {
      stroke(173, 255, 47, 80);
      fill(173, 255, 47, 80);
    }
    if ((strand1.equals("-")) && (strand2.equals("-"))) {
      quad(x2_symb1+genome_height+offset1, y1+(genome_height/4), x1_symb1+genome_height+offset1, y1+(genome_height/4), x1_symb2+genome_height+offset2, y2+(genome_height/4), x2_symb2+genome_height+offset2, y2+(genome_height/4));
    }

    if ((strand1.equals("+")) && (strand2.equals("+"))) {
      quad(x2_symb1+genome_height+offset1, y1-(genome_height/4), x1_symb1+genome_height+offset1, y1-(genome_height/4), x1_symb2+genome_height+offset2, y2-(genome_height/4), x2_symb2+genome_height+offset2, y2-(genome_height/4));
    }

    if ((strand1.equals("+")) && (strand2.equals("-"))) {
      quad(x2_symb1+genome_height+offset1, y1-(genome_height/4), x1_symb1+genome_height+offset1, y1-(genome_height/4), x1_symb2+genome_height+offset2, y2+(genome_height/4), x2_symb2+genome_height+offset2, y2+(genome_height/4));
    }
    if ((strand1.equals("-")) && (strand2.equals("+"))) {
      quad(x2_symb1+genome_height+offset1, y1+(genome_height/4), x1_symb1+genome_height+offset1, y1+(genome_height/4), x1_symb2+genome_height+offset2, y2-(genome_height/4), x2_symb2+genome_height+offset2, y2-(genome_height/4));
    }
  }
}

void draw() {

  drawConnections(connections[0], genome[0], genome[1], 2*genome_height, 4*genome_height); //csv file, genome size 1, genome size 2, y1, y2
  drawConnections(connections[1], genome[1], genome[2], 4*genome_height, 6*genome_height);
  drawConnections(connections[2], genome[2], genome[3], 6*genome_height, 8*genome_height);
  drawConnections(connections[3], genome[3], genome[4], 8*genome_height, 10*genome_height);

  drawGenes(genes[0], genome[0], 2*genome_height);
  drawGenes(genes[1], genome[1], 4*genome_height);
  drawGenes(genes[2], genome[2], 6*genome_height);
  drawGenes(genes[3], genome[3], 8*genome_height);
  drawGenes(genes[4], genome[4], 10*genome_height);

  //TODO: for loop in the draw function

  save("genome_alignment.tif");
  endRecord();
}