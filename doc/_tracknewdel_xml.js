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
'deprecated.xml',
'tracknewdel.xml'
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
'include_deprecated.xml',
'fundeprecated.xml',
'omp_max_thread.xml',
'tracknewdel.xml',
'omp_alloc.xml'
];
var list_down0 = [
'tracknewdel.cpp.xml'
];
var list_current0 = [
'tracknewdel.xml#Deprecated',
'tracknewdel.xml#Syntax',
'tracknewdel.xml#Purpose',
'tracknewdel.xml#Include',
'tracknewdel.xml#file',
'tracknewdel.xml#line',
'tracknewdel.xml#oldptr',
'tracknewdel.xml#newlen',
'tracknewdel.xml#head newptr',
'tracknewdel.xml#ncopy',
'tracknewdel.xml#TrackNewVec',
'tracknewdel.xml#TrackNewVec.Macro',
'tracknewdel.xml#TrackNewVec.Previously Deprecated',
'tracknewdel.xml#TrackDelVec',
'tracknewdel.xml#TrackDelVec.Macro',
'tracknewdel.xml#TrackDelVec.Previously Deprecated',
'tracknewdel.xml#TrackExtend',
'tracknewdel.xml#TrackExtend.Macro',
'tracknewdel.xml#TrackExtend.Previously Deprecated',
'tracknewdel.xml#TrackCount',
'tracknewdel.xml#TrackCount.Macro',
'tracknewdel.xml#TrackCount.Previously Deprecated',
'tracknewdel.xml#Multi-Threading',
'tracknewdel.xml#Example'
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
