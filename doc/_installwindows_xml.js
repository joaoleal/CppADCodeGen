var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'install.xml',
'installwindows.xml'
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
'installunix.xml',
'pkgconfig.xml',
'installwindows.xml'
];
var list_current0 = [
'installwindows.xml#Cygwin',
'installwindows.xml#Download',
'installwindows.xml#Download.Subversion',
'installwindows.xml#Download.Web Link',
'installwindows.xml#Download.Unix Tar Files',
'installwindows.xml#Getting Started',
'installwindows.xml#Introduction',
'installwindows.xml#Examples and Testing',
'installwindows.xml#More Correctness Testing',
'installwindows.xml#Printing During Forward Mode',
'installwindows.xml#CppAD Speed Test',
'installwindows.xml#Double Speed Test',
'installwindows.xml#Speed Utility Example'
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
