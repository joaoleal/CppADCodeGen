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
'exp_eps.xml'
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
'get_started.cpp.xml',
'exp_2.xml',
'exp_eps.xml',
'exp_apx_main.cpp.xml'
];
var list_down0 = [
'exp_eps.hpp.xml',
'exp_eps.cpp.xml',
'exp_eps_for0.xml',
'exp_eps_for1.xml',
'exp_eps_rev1.xml',
'exp_eps_for2.xml',
'exp_eps_rev2.xml',
'exp_eps_cppad.xml'
];
var list_current0 = [
'exp_eps.xml#Syntax',
'exp_eps.xml#Purpose',
'exp_eps.xml#Mathematical Function',
'exp_eps.xml#include',
'exp_eps.xml#x',
'exp_eps.xml#epsilon',
'exp_eps.xml#y',
'exp_eps.xml#Type',
'exp_eps.xml#Implementation',
'exp_eps.xml#Test',
'exp_eps.xml#Exercises'
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
