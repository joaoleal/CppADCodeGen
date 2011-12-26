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
'exp_eps.xml',
'exp_eps_for1.xml'
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
'exp_eps.hpp.xml',
'exp_eps.cpp.xml',
'exp_eps_for0.xml',
'exp_eps_for1.xml',
'exp_eps_rev1.xml',
'exp_eps_for2.xml',
'exp_eps_rev2.xml',
'exp_eps_cppad.xml'
];
var list_down0 = [
'exp_eps_for1.cpp.xml'
];
var list_current0 = [
'exp_eps_for1.xml#First Order Expansion',
'exp_eps_for1.xml#Mathematical Form',
'exp_eps_for1.xml#Operation Sequence',
'exp_eps_for1.xml#Operation Sequence.Index',
'exp_eps_for1.xml#Operation Sequence.Operation',
'exp_eps_for1.xml#Operation Sequence.Zero Order',
'exp_eps_for1.xml#Operation Sequence.Derivative',
'exp_eps_for1.xml#Operation Sequence.First Order',
'exp_eps_for1.xml#Operation Sequence.Sweep',
'exp_eps_for1.xml#Return Value',
'exp_eps_for1.xml#Verification',
'exp_eps_for1.xml#Exercises'
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
