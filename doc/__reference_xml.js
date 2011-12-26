var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'_reference.xml'
];
var list_down1 = [
'install.xml',
'introduction.xml',
'ad.xml',
'adfun.xml',
'multi_thread.xml',
'library.xml',
'cppad_ipopt_nlp.xml',
'example.xml',
'preprocessor.xml',
'appendix.xml'
];
var list_current0 = [
'_reference.xml#A',
'_reference.xml#B',
'_reference.xml#C',
'_reference.xml#D',
'_reference.xml#E',
'_reference.xml#F',
'_reference.xml#G',
'_reference.xml#H',
'_reference.xml#I',
'_reference.xml#J',
'_reference.xml#L',
'_reference.xml#M',
'_reference.xml#N',
'_reference.xml#O',
'_reference.xml#P',
'_reference.xml#R',
'_reference.xml#S',
'_reference.xml#T',
'_reference.xml#U',
'_reference.xml#V',
'_reference.xml#W'
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
