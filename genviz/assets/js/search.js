var formatSearchTerm = (input_gene, input_orgn) => {
	// TODO: Validate gene input
	var gene = input_gene.val()
	var organism = (input_orgn ? input_orgn.val() : null) || 'Homo sapiens'
	
	return `"${gene}"[gene] AND "${organism}"[orgn]`
}