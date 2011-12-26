var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'introduction.xml',
'exp_2.xml',
'exp_2_for0.xml'
];
var list_down3 = [
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
var list_down2 = [
'get_started.cpp.xml',
'exp_2.xml',
'exp_eps.xml',
'exp_apx_main.cpp.xml'
];
var list_down1 = [
'exp_2.hpp.xml',
'exp_2.cpp.xml',
'exp_2_for0.xml',
'exp_2_for1.xml',
'exp_2_rev1.xml',
'exp_2_for2.xml',
'exp_2_rev2.xml',
'exp_2_cppad.xml'
];
var list_down0 = [
'exp_2_for0.cpp.xml'
];
var list_current0 = [
'exp_2_for0.xml#Mathematical Form',
'exp_2_for0.xml#Zero Order Expansion',
'exp_2_for0.xml#Operation Sequence',
'exp_2_for0.xml#Operation Sequence.Index',
'exp_2_for0.xml#Operation Sequence.Code',
'exp_2_for0.xml#Operation Sequence.Operation',
'exp_2_for0.xml#Operation Sequence.Zero Order',
'exp_2_for0.xml#Operation Sequence.Sweep',
'exp_2_for0.xml#Return Value',
'exp_2_for0.xml#Verification',
'exp_2_for0.xml#Exercises'
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
