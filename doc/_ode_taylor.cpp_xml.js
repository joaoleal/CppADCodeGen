var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'example.xml',
'general.xml',
'ode_taylor.cpp.xml'
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
'general.xml',
'exampleutility.xml',
'listallexamples.xml',
'test_vector.xml'
];
var list_down1 = [
'ad_fun.cpp.xml',
'ad_in_c.cpp.xml',
'conj_grad.cpp.xml',
'hesminordet.cpp.xml',
'hesludet.cpp.xml',
'interface2c.cpp.xml',
'jacminordet.cpp.xml',
'jacludet.cpp.xml',
'mul_level.xml',
'odestiff.cpp.xml',
'ode_taylor.cpp.xml',
'ode_taylor_adolc.cpp.xml',
'stackmachine.cpp.xml'
];
var list_current0 = [
'ode_taylor.cpp.xml#Purpose',
'ode_taylor.cpp.xml#ODE',
'ode_taylor.cpp.xml#ODE Solution',
'ode_taylor.cpp.xml#Derivative of ODE Solution',
'ode_taylor.cpp.xml#Taylor\'s Method Using AD'
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
