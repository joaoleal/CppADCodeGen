var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'appendix.xml',
'speed.xml',
'speed_utility.xml',
'ode_evaluate.xml'
];
var list_down3 = [
'faq.xml',
'speed.xml',
'theory.xml',
'glossary.xml',
'bib.xml',
'bugs.xml',
'wishlist.xml',
'whats_new.xml',
'deprecated.xml',
'license.xml'
];
var list_down2 = [
'speed_main.xml',
'speed_utility.xml',
'speed_double.xml',
'speed_adolc.xml',
'speed_cppad.xml',
'speed_fadbad.xml',
'speed_sacado.xml'
];
var list_down1 = [
'uniform_01.xml',
'det_of_minor.xml',
'det_by_minor.xml',
'det_by_lu.xml',
'det_33.xml',
'det_grad_33.xml',
'mat_sum_sq.xml',
'ode_evaluate.xml',
'sparse_evaluate.xml'
];
var list_down0 = [
'ode_evaluate.cpp.xml',
'ode_evaluate.hpp.xml'
];
var list_current0 = [
'ode_evaluate.xml#Syntax',
'ode_evaluate.xml#Purpose',
'ode_evaluate.xml#Inclusion',
'ode_evaluate.xml#Float',
'ode_evaluate.xml#x',
'ode_evaluate.xml#m',
'ode_evaluate.xml#m = 0',
'ode_evaluate.xml#m = 1',
'ode_evaluate.xml#fm',
'ode_evaluate.xml#fm.Function',
'ode_evaluate.xml#fm.Gradient',
'ode_evaluate.xml#Example',
'ode_evaluate.xml#Source Code'
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
