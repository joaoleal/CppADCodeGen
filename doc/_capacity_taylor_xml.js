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
'forward.xml',
'capacity_taylor.xml'
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
'forwardzero.xml',
'forwardone.xml',
'forwardany.xml',
'size_taylor.xml',
'comparechange.xml',
'capacity_taylor.xml',
'forward.cpp.xml'
];
var list_current0 = [
'capacity_taylor.xml#Syntax',
'capacity_taylor.xml#Syntax.See Also',
'capacity_taylor.xml#Purpose',
'capacity_taylor.xml#f',
'capacity_taylor.xml#c',
'capacity_taylor.xml#Discussion',
'capacity_taylor.xml#Discussion.Pre-Allocating Memory',
'capacity_taylor.xml#Discussion.Freeing Memory',
'capacity_taylor.xml#Discussion.Original State',
'capacity_taylor.xml#Example'
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
