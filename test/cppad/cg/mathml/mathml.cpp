/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2015 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

#include <iostream>
#include <fstream>

#include <cppad/cg/cppadcg.hpp>
#include <cppad/cg/mathml/mathml.hpp>
#include <gtest/gtest.h>

using namespace CppAD;
using namespace CppAD::cg;

TEST(CppADCGLatexTest, latex) {
    // use a special object for source code generation
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;

    // independent variable vector
    CppAD::vector<ADCG> x(2);
    x[0] = 2.;
    x[1] = 3.;
    Independent(x);

    // dependent variable vector 
    CppAD::vector<ADCG> y(8);

    // the model
    ADCG a = x[0] / 1. + x[1] * x[1];
    ADCG b = a / 2e-6;
    y[0] = b + 1 / (sign(b)*5 * a);
    y[1] = x[1];
    y[2] = CondExpLt(ADCG(1.0), x[0], x[1], b);
    y[3] = CondExpLe(x[0], ADCG(2.0), x[1], b);
    y[4] = CondExpEq(x[0], x[1], x[1], b);
    ADCG c = CondExpGe(ADCG(3.0), x[0], a, b);
    y[5] = CondExpGt(ADCG(4.0), x[0], ADCG(5.0), c);
    y[6] = 5 * pow(4, x[0]);
    y[7] = 3;

    ADFun<CGD> fun(x, y); // the model tape

    /**
     * start the special steps for source code generation
     */
    CodeHandler<double> handler;

    CppAD::vector<CGD> indVars(2);
    handler.makeVariables(indVars);

    CppAD::vector<CGD> vals = fun.Forward(0, indVars);

    LanguageMathML<double> langMathML;
    LangMathMLDefaultVariableNameGenerator<double> nameGen;

    langMathML.setSaveVariableRelations(true);

    // add some additional code to select variables
    langMathML.setStyle(langMathML.getStyle() + "\n.selected{background-color: #ccc;}"
                                                "\n.faded{\n"
                                                "    opacity: 0.3;\n"
                                                "    filter: alpha(opacity=30); /* For IE8 and earlier */\n"
                                                "}");

    // use block display
    langMathML.setEquationMarkup("<math display=\"block\" class=\"equation\">", "</math>");

    // use inline display
    //langMathML.setEquationMarkup("<math>", "</math><br/>");

    // use MathJax (and align to the left)
    langMathML.setHeadExtraMarkup("<script type=\"text/x-mathjax-config\">\n"
                                  //"MathJax.Hub.Config({    MMLorHTML: { prefer: { Firefox: \"MML\" } }  });\n" // use this to define a prefered browser renderer
                                  "MathJax.Hub.Config({\n"
                                  "    jax: [\"input/TeX\",\"output/HTML-CSS\"],\n"
                                  "    displayAlign: \"left\"\n"
                                  "});\n"
                                  "</script>\n"
                                  "<script type=\"text/javascript\" src=\"https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script>");

    langMathML.setJavascript("var selected = [];\n"
                             "\n"
                             "function contains(arr, o) {\n"
                             "    if(usages === null || usages === undefined) {\n"
                             "        return false;\n"
                             "    }\n"
                             "    var l = arr.length;\n"
                             "    for (var i = 0; i < l; i++) {\n"
                             "        if (arr[i] === o) {\n"
                             "            return true;\n"
                             "        }\n"
                             "    }\n"
                             "    return false;\n"
                             "}\n"
                             "\n"
                             "function findEquation(el) {\n"
                             "    while (el != document) {\n"
                             "        if(el.classList.contains('equation')) {\n"
                             "            return el;\n"
                             "        }\n"
                             "        el = el.parentNode;\n"
                             "    }\n"
                             "    return null;\n"
                             "}\n"
                             "\n"
                             "function isBranch(eq) {\n"
                             "    var el = eq.parentNode;\n"
                             "    while (el != document && el.id != 'algorithm') {\n"
                             "        if(el.classList.contains('condBody')) {\n"
                             "            return true;\n"
                             "        }\n"
                             "        el = el.parentNode;\n"
                             "    }\n"
                             "    return false;\n"
                             "}\n"
                             "\n"
                             "function showEquations(eqId) {\n"
                             "    var deps = var2dep[eqId];\n"
                             "    if(deps === undefined || deps === null) {\n"
                             "        return;\n"
                             "    }\n"
                             "    l = deps.length;\n"
                             "    for(var i = 0; i < l; i++) {\n"
                             "        var id = deps[i];\n"
                             "        var el = document.getElementById(\"v\" + id);\n"
                             "        if(el !== null) {\n"
                             "            var eq = findEquation(el);\n"
                             "            if(eq !== null && eq.classList.contains('faded')) {    \n"
                             "                eq.classList.remove('faded');\n"
                             "                showEquations(id);\n"
                             "            }\n"
                             "        }\n"
                             "    }\n"
                             "}\n"
                             "\n"
                             "function hideEquationForIds(ids, visibleId) {\n"
                             "    for(var co in ids) {\n"
                             "        if(visibleId === ids[co])\n"
                             "            continue;\n"
                             "        var el = document.getElementById(ids[co]);\n"
                             "        if(el !== null && el !== undefined) {\n"
                             "            var eq = findEquation(el);\n"
                             "            if(eq !== null && eq !== undefined)\n"
                             "                eq.classList.add('faded');\n"
                             "        }\n"
                             "    }\n"
                             "}\n"
                             "function clearAllClass(className) {\n"
                             "    var list = document.getElementsByClassName(className);\n"
                             "    if(list !== undefined) {\n"
                             "        while(list.length > 0){\n"
                             "            list[0].classList.remove(className);\n"
                             "        }\n"
                             "    }\n"
                             "}"
                             "\n"
                             "function clickHandler(e) {\n"
                             "    var t = e.target;\n"
                             "    \n"
                             "    clearAllClass('selected');\n"
                             "    clearAllClass('faded');\n"
                             "   \n"
                             "    while (t != document) {\n"
                             "        if (t.id !== null && t.id !== \"\" && t.id.charAt(0) == 'v') {\n"
                             "            var baseId = t.id.split('_')[0];\n"
                             "            var idval = baseId.substring(1);\n"
                             "            var el = document.getElementById(baseId);\n"
                             "            var n = 0;\n"
                             "            while (el !== null) {\n"
                             "                el.classList.add('selected');\n"
                             "                n++;\n"
                             "                el = document.getElementById(baseId + '_' + n);\n"
                             "            }\n"
                             "            \n"
                             "            // fade other equations which do not use this variable\n"
                             "            usages = dep2var[idval];\n"
                             "            for (var i in var2dep) {\n"
                             "                if(i === idval)\n"
                             "                    continue;\n"
                             "                var vi = parseInt(i);\n"
                             "                if(!contains(usages, vi)) {\n"
                             "                    var el2 = document.getElementById(\"v\" + i);\n"
                             "                    n = 0;\n"
                             "                    while(el2 !== null) {\n"
                             "                        var eq = findEquation(el2);\n"
                             "                        if(eq === null)\n"
                             "                            break;\n"
                             "                        eq.classList.add('faded');\n"
                             "                        if(!isBranch(eq))\n"
                             "                            break;\n"
                             "                        n++;\n"
                             "                        el2 = document.getElementById(\"v\" + i + '_' + n);\n"
                             "                    }\n"
                             "                }\n"
                             "            }\n"
                             "            \n"
                             "            // unfade equations used to create that variable\n"
                             "            showEquations(idval);\n"
                             "            \n"
                             "            hideEquationForIds(depConst, t.id);\n"
                             "            hideEquationForIds(depIsVar, t.id);\n"
                             "            break;\n"
                             "        }\n"
                             "        t = t.parentNode;\n"
                             "    }\n"
                             "}\n"
                             "\n"
                             "document.addEventListener('DOMContentLoaded', function () {\n"
                             "    document.getElementById('algorithm').onclick = clickHandler;\n"
                             "}, false);");

    // create the HMTL file
    std::ofstream htmlFile;
    htmlFile.open("algorithm.html");

    handler.setReuseVariableIDs(false);
    handler.generateCode(htmlFile, langMathML, vals, nameGen);

    htmlFile.close();
}
