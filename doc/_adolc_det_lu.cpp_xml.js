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
'speed_adolc.xml',
'adolc_det_lu.cpp.xml'
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
'adolc_det_minor.cpp.xml',
'adolc_det_lu.cpp.xml',
'adolc_mat_mul.cpp.xml',
'adolc_ode.cpp.xml',
'adolc_poly.cpp.xml',
'adolc_sparse_hessian.cpp.xml',
'adolc_sparse_jacobian.cpp.xml'
];
var list_current0 = [
'adolc_det_lu.cpp.xml#Specifications',
'adolc_det_lu.cpp.xml#Implementation'
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
