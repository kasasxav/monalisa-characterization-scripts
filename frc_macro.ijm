function action(input, output, filename) {
	open(input + filename);
	selectWindow(filename);
	run("Stack to Images");
	run("FRC Calculation...", "image_1=04c-0001 image_2=04c-0002 resolution=[Fixed 1/7] display");
	Plot.getValues(x, y);
	for (i=0; i<x.length; i++) {
		setResult("X", i, x[i]);
		setResult("Y", i, y[i]);
	}
	saveAs("Results", output+filename+".txt");
	run("Close All");
}

input = "C:/Users/xavie/Documents/conf/";
output = "C:/Users/xavie/Documents/conf_output/";
list = getFileList(input);
for (i = 0; i < list.length; i++){
        action(input, output, list[i]);
}