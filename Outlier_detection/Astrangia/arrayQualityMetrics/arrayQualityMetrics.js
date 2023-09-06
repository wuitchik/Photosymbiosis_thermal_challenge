// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false ];
var arrayMetadata    = [ [ "1", "AA2", "1.80838514054066", "heat_1", "White", "Heat", "3", "A" ], [ "2", "AA3", "1.1937693740877", "cold_1", "White", "Cold", "1", "A" ], [ "3", "AA4", "1.04777719237146", "control_1", "White", "Control", "5", "A" ], [ "4", "AB1", "1.83209680834956", "control_2", "Brown", "Control", "5", "B" ], [ "5", "AB2", "0.958808748100575", "cold_1", "Brown", "Cold", "1", "B" ], [ "6", "AB3", "0.422664331708805", "heat_1", "Brown", "Heat", "3", "B" ], [ "7", "AC1", "2.2749291969112", "heat_2", "White", "Heat", "5", "C" ], [ "8", "AC3", "2.8576230591922", "control_2", "White", "Control", "5", "C" ], [ "9", "AD4", "1.49728064087956", "control_3", "White", "Control", "5", "D" ], [ "10", "AD5", "0.43603777952792", "heat_2", "White", "Heat", "5", "D" ], [ "11", "AD6", "1.79739738018844", "cold_1", "White", "Cold", "1", "D" ], [ "12", "AE2", "1.42946449921428", "cold_3", "White", "Cold", "1", "E" ], [ "13", "AE5", "0.893635681986263", "heat_1", "White", "Heat", "5", "E" ], [ "14", "AE6", "1.36617479428685", "control_2", "White", "Control", "3", "E" ], [ "15", "AF2", "0.0146246657103813", "control_2", "Brown", "Control", "5", "F" ], [ "16", "AF3", "0.998964894983688", "heat_3", "Brown", "Heat", "5", "F" ], [ "17", "AF5", "1.88529772955377", "cold_1", "Brown", "Cold", "1", "F" ], [ "18", "AG2", "2.06485044784716", "control_3", "White", "Control", "5", "G" ], [ "19", "AG3", "2.35553046926386", "cold_3", "White", "Cold", "1", "G" ], [ "20", "AH1", "1.16014659235784", "heat_3", "White", "Heat", "5", "H" ], [ "21", "AH2", "0.526635560992864", "control_1", "White", "Control", "5", "H" ], [ "22", "AH3", "1.13157008442994", "cold_1", "White", "Cold", "1", "H" ], [ "23", "AI1", "1.90269262031657", "control_3", "Brown", "Control", "5", "I" ], [ "24", "AI2", "0.576862123088409", "heat_3", "Brown", "Heat", "2", "I" ], [ "25", "AI5", "0.99436165528483", "cold_2", "Brown", "Cold", "2", "I" ], [ "26", "AJ2", "1.21640948091622", "cold_3", "Brown", "Cold", "1", "J" ], [ "27", "AJ3", "0.908415194561251", "heat_2", "Brown", "Heat", "5", "J" ], [ "28", "AJ4", "1.1056284592401", "control_1", "Brown", "Control", "5", "J" ], [ "29", "AK3", "0.849598008723261", "cold_2", "White", "Cold", "1", "K" ], [ "30", "AK5", "1.56401736164768", "control_2", "White", "Control", "5", "K" ], [ "31", "AL1", "1.17129123202391", "cold_2", "Brown", "Cold", "1", "L" ], [ "32", "AL2", "0.570002073792678", "heat_2", "Brown", "Heat", "4", "L" ], [ "33", "AM1", "0.591065203922657", "cold_1", "Brown", "Cold", "1", "M" ], [ "34", "AM2", "1.45470324301181", "heat_1", "Brown", "Heat", "4", "M" ], [ "35", "AM3", "1.64510366774562", "control_1", "Brown", "Control", "5", "M" ], [ "36", "AN1", "1.37230679371832", "heat_3", "Brown", "Heat", "1", "N" ], [ "37", "AN2", "2.01486380521999", "cold_3", "Brown", "Cold", "1", "N" ], [ "38", "AN3", "0.823753444727633", "control_3", "Brown", "Control", "5", "N" ], [ "39", "AP1", "1.63219504613154", "cold_3", "White", "Cold", "1", "P" ], [ "40", "AP2", "0.932056629738341", "control_1", "White", "Control", "5", "P" ], [ "41", "AP4", "1.11662328848137", "heat_3", "White", "Heat", "4", "P" ], [ "42", "AS3", "0.233879167576276", "heat_1", "Brown", "Heat", "2", "S" ], [ "43", "AS5", "0.923025490832624", "cold_2", "Brown", "Cold", "1", "S" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
