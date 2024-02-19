const fs = require('fs');
const genbankParser = require('genbank-parser');
 
const genbank = fs.readFileSync('./GCA_003054575.1.gb', 'utf-8');
const result = genbankParser(genbank);