var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'cppad_ipopt_nlp.xml'
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
var list_down0 = [
'cppad_ipopt_windows.xml',
'ipopt_get_started.cpp.xml',
'cppad_ipopt_ode.xml',
'ipopt_ode_speed.cpp.xml'
];
var list_current0 = [
'cppad_ipopt_nlp.xml#Syntax',
'cppad_ipopt_nlp.xml#Purpose',
'cppad_ipopt_nlp.xml#ipopt_library_paths',
'cppad_ipopt_nlp.xml#fg(x)',
'cppad_ipopt_nlp.xml#fg(x).Index Vector',
'cppad_ipopt_nlp.xml#fg(x).Projection',
'cppad_ipopt_nlp.xml#fg(x).Injection',
'cppad_ipopt_nlp.xml#fg(x).Representation',
'cppad_ipopt_nlp.xml#Simple Representation',
'cppad_ipopt_nlp.xml#SizeVector',
'cppad_ipopt_nlp.xml#NumberVector',
'cppad_ipopt_nlp.xml#ADNumber',
'cppad_ipopt_nlp.xml#ADVector',
'cppad_ipopt_nlp.xml#n',
'cppad_ipopt_nlp.xml#m',
'cppad_ipopt_nlp.xml#x_i',
'cppad_ipopt_nlp.xml#x_l',
'cppad_ipopt_nlp.xml#x_u',
'cppad_ipopt_nlp.xml#g_l',
'cppad_ipopt_nlp.xml#g_u',
'cppad_ipopt_nlp.xml#fg_info',
'cppad_ipopt_nlp.xml#fg_info.fg_info.number_functions',
'cppad_ipopt_nlp.xml#fg_info.fg_info.eval_r',
'cppad_ipopt_nlp.xml#fg_info.fg_info.retape',
'cppad_ipopt_nlp.xml#fg_info.fg_info.domain_size',
'cppad_ipopt_nlp.xml#fg_info.fg_info.range_size',
'cppad_ipopt_nlp.xml#fg_info.fg_info.number_terms',
'cppad_ipopt_nlp.xml#fg_info.fg_info.index',
'cppad_ipopt_nlp.xml#solution',
'cppad_ipopt_nlp.xml#solution.status',
'cppad_ipopt_nlp.xml#solution.x',
'cppad_ipopt_nlp.xml#solution.z_l',
'cppad_ipopt_nlp.xml#solution.z_u',
'cppad_ipopt_nlp.xml#solution.g',
'cppad_ipopt_nlp.xml#solution.lambda',
'cppad_ipopt_nlp.xml#solution.obj_value',
'cppad_ipopt_nlp.xml#Visual Studio',
'cppad_ipopt_nlp.xml#Example',
'cppad_ipopt_nlp.xml#Wish List'
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
