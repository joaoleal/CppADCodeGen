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
'seq_property.xml'
];
var list_down2 = [
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
var list_down1 = [
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
var list_down0 = [
'seq_property.cpp.xml'
];
var list_current0 = [
'seq_property.xml#Syntax',
'seq_property.xml#Syntax.See Also',
'seq_property.xml#Purpose',
'seq_property.xml#f',
'seq_property.xml#Domain',
'seq_property.xml#Range',
'seq_property.xml#Parameter',
'seq_property.xml#size_var',
'seq_property.xml#size_par',
'seq_property.xml#size_op',
'seq_property.xml#size_op_arg',
'seq_property.xml#size_text',
'seq_property.xml#size_VecAD',
'seq_property.xml#size_op_seq',
'seq_property.xml#Example'
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
