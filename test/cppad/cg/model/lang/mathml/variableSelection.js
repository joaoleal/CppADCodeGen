function contains(arr, o) {
    if (usages === null || usages === undefined) {
        return false;
    }
    var l = arr.length;
    for (var i = 0; i < l; i++) {
        if (arr[i] === o) {
            return true;
        }
    }
    return false;
}

function findEquation(el) {
    if (el.classList.contains('indep'))
        return null;
    while (el != document) {
        if (el.classList.contains('equation')) {
            return el;
        }
        el = el.parentNode;
    }
    return null;
}

function isBranch(eq) {
    var el = eq.parentNode;
    while (el != document && el.id != 'algorithm') {
        if (el.classList.contains('condBody')) {
            return true;
        }
        el = el.parentNode;
    }
    return false;
}

function showEquations(eqId, level) {
    var el = document.getElementById('v' + eqId);
    if (el !== null) {
        var eq = findEquation(el);
        if (eq !== null && eq !== undefined) {
            eq.classList.remove('faded');
            if (!eq.classList.contains('depEq')) {
                eq.classList.add('depEq');
                if (level > 1)
                    eq.classList.add('faded2');
            }
        }
    }

    var deps = var2dep[eqId];
    if (deps === undefined || deps === null) {
        return;
    }
    for (var i = 0; i < deps.length; i++) {
        var id = deps[i];
        showEquations(id, level + 1);
    }
}

function hideEquationForIds(ids, visibleId) {
    for (var co in ids) {
        if (visibleId === ids[co])
            continue;
        var el = document.getElementById(ids[co]);
        if (el !== null && el !== undefined) {
            var eq = findEquation(el);
            if (eq !== null && eq !== undefined)
                eq.classList.add('faded');
        }
    }
}
function clearAllClass(className) {
    var list = document.getElementsByClassName(className);
    if (list !== undefined) {
        while (list.length > 0) {
            list[0].classList.remove(className);
        }
    }
}

function clickHandler(e) {
    var t = e.target;

    clearAllClass('selectedProp');
    clearAllClass('faded');
    clearAllClass('faded2');
    clearAllClass('depEq');

    while (t != document) {
        if (t.id !== null && t.id !== '' && t.id.charAt(0) == 'v') {
            var baseId = t.id.split('_')[0];
            var idval = baseId.substring(1);
            var el = document.getElementById(baseId);
            var n = 0;
            while (el !== null) {
                el.classList.add('selectedProp');
                n++;
                el = document.getElementById(baseId + '_' + n);
            }

            // fade other equations which do not use this variable
            usages = dep2var[idval];
            for (var i in var2dep) {
                if (i === idval)
                    continue;
                var vi = parseInt(i);
                if (!contains(usages, vi)) {
                    var el2 = document.getElementById('v' + i);
                    n = 0;
                    while (el2 !== null) {
                        var eq = findEquation(el2);
                        if (eq === null)
                            break;
                        eq.classList.add('faded');
                        if (!isBranch(eq))
                            break;
                        n++;
                        el2 = document.getElementById('v' + i + '_' + n);
                    }
                }
            }

            // unfade equations used to create that variable
            showEquations(idval, 0);

            hideEquationForIds(depConst, t.id);
            hideEquationForIds(depIsVar, t.id);
            break;
        }
        t = t.parentNode;
    }
}

document.addEventListener('DOMContentLoaded', function () {
    document.getElementById('algorithm').onclick = clickHandler;
}, false);