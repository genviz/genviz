Handlebars.registerHelper("jsonStringify", function(value, options) {
    return JSON.stringify(value);
});

Handlebars.registerHelper("last", function(arr, options) {
    return arr[arr.length - 1];
});

Handlebars.registerHelper("inc", function(value, options) {
    return parseInt(value) + 1;
});

Handlebars.registerHelper("dec", function(value, options) {
    return parseInt(value) - 1;
});

Handlebars.registerHelper("sum", function(a, b, options) {
    return parseInt(a) + parseInt(b);
});
Handlebars.registerHelper("substract", function(a, b, options) {
    return parseInt(a) - parseInt(b);
});

Handlebars.registerHelper("multiply", function(a, b, options) {
	return parseInt(a) * parseInt(b);
});

Handlebars.registerHelper("eq", function(a, b, options) {
	return a == b;
});

Handlebars.registerHelper("divisible", function(a, b, options) {
    return a % b == 0;
});

Handlebars.registerHelper("mod", function(a, b, options) {
    return a % b;
});

Handlebars.registerHelper('times', function(n, block) {
    var accum = '';
    for(var i = 0; i < n; ++i)
        accum += block.fn(i);
    return accum;
});

Handlebars.registerHelper('times1', function(n, block) {
    var accum = '';
    for(var i = 1; i <= n; ++i)
        accum += block.fn(i);
    return accum;
});

Handlebars.registerHelper('eachChar', function(str, block) {
    var accum = ''
    for(var i = 0; i < str.length; ++i)
        accum += block.fn(str[i]);
    return accum;
});


Handlebars.registerHelper("sourceOffset", function(vars, row_i, basesPerRow, offset, options) {
    console.log("Vars", vars, (row_i+1)*basesPerRow, vars[vars.length - 1].end, offset, vars[vars.length - 1].end - offset)
    return (row_i+1) * basesPerRow - (vars[vars.length - 1].end - offset) + 1;
});
