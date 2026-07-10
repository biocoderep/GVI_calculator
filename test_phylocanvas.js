const Phylocanvas = require('./phylocanvas.js');
console.log("Keys in Phylocanvas:", Object.keys(Phylocanvas));
console.log("createTree typeof:", typeof Phylocanvas.createTree);
if (Phylocanvas.default) {
    console.log("Keys in Phylocanvas.default:", Object.keys(Phylocanvas.default));
}
