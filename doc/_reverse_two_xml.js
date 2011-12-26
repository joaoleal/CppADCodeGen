var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'adfun.xml',
'funeval.xml',
'reverse.xml',
'reverse_two.xml'
];
var list_down3 = [
'independent.xml',
'funconstruct.xml',
'dependent.xml',
'abort_recording.xml',
'seq_property.xml',
'funeval.xml',
'drivers.xml',
'funcheck.xml',
'optimize.xml'
];
var list_down2 = [
'forward.xml',
'reverse.xml',
'sparse.xml'
];
var list_down1 = [
'reverse_one.xml',
'reverse_two.xml',
'reverse_any.xml'
];
var list_down0 = [
'reverse_two.cpp.xml',
'hestimesdir.cpp.xml'
];
var list_current0 = [
'reverse_two.xml#Syntax',
'reverse_two.xml#Purpose',
'reverse_two.xml#x^(k)',
'reverse_two.xml#W',
'reverse_two.xml#f',
'reverse_two.xml#w',
'reverse_two.xml#dw',
'reverse_two.xml#dw.First Order Partials',
'reverse_two.xml#dw.Second Order Partials',
'reverse_two.xml#Vector',
'reverse_two.xml#Hessian Times Direction',
'reverse_two.xml#Example'
];
function choose_across0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_across0[index-1];
}
function choose_up0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_up0[index-1];
}
function choose_down3(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down3[index-1];
}
function choose_down2(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down2[index-1];
}
function choose_down1(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down1[index-1];
}
function choose_down0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down0[index-1];
}
function choose_current0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_current0[index-1];
}
