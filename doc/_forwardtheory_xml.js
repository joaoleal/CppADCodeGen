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
'theory.xml',
'forwardtheory.xml'
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
var list_down1 = [
'forwardtheory.xml',
'reversetheory.xml',
'reverse_identity.xml'
];
var list_down0 = [
'expforward.xml',
'logforward.xml',
'sqrtforward.xml',
'sincosforward.xml',
'atanforward.xml',
'asinforward.xml',
'acosforward.xml',
'tan_forward.xml'
];
var list_current0 = [
'forwardtheory.xml#Taylor Notation',
'forwardtheory.xml#Binary Operators',
'forwardtheory.xml#Binary Operators.Addition',
'forwardtheory.xml#Binary Operators.Subtraction',
'forwardtheory.xml#Binary Operators.Multiplication',
'forwardtheory.xml#Binary Operators.Division',
'forwardtheory.xml#Standard Math Functions',
'forwardtheory.xml#Standard Math Functions.Differential Equation',
'forwardtheory.xml#Standard Math Functions.Taylor Coefficients Recursion Formula',
'forwardtheory.xml#Standard Math Functions.Cases that Apply Recursion Above',
'forwardtheory.xml#Standard Math Functions.Special Cases'
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
