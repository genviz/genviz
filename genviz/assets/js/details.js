function scrollToAnchor(id){
    var tag = $("#"+id);
    $('html,body').animate({scrollTop: tag.offset().top});
}

function checkInView(elem, container)
{
    var contHeight = container.height()
    var contTop = container.scrollTop()
    var contBottom = contTop + contHeight
    var elemTop = $(elem).offset().top - container.offset().top
    var elemBottom = elemTop + $(elem).height()
    var isTotal = (elemTop >= 0 && elemBottom <= contHeight)
    // var isPart = ((elemTop < 0 && elemBottom > 0 ) || (elemTop > 0 && elemTop <= container.height()));

    return  isTotal
}

function getFormData($form){
    var unindexed_array = $form.serializeArray();
    var indexed_array = {};

    $.map(unindexed_array, function(n, i){
        indexed_array[n['name']] = n['value'];
    });

    return indexed_array;
}

function plotGeneFeatures(features) {
	// Dummy timestamp to use non-date values in timeline vis.js chart
	var dummyTs = +new Date()

	// Build vis.js dataset
	var minLocation = Infinity;
	var maxLocation = -Infinity;
	var exons_dp = []
	var prev_end = -1
	var row = 0
	Object.keys(features).forEach(function(feature_type, _) {
		features[feature_type].forEach(function(feature, i) {
			start = feature.location[0]
			end = feature.location[1]
			minLocation = Math.min(minLocation, start)
			maxLocation = Math.max(maxLocation, end)
			// Start and end date are the locations of the features + a dummy ts
			// so it thinks it's a date and can be used with timeline
			exons_dp.push({
				start: start + dummyTs,
				end: end + dummyTs,
				content: `${feature_type} ${i+1}`, group: feature_type,
				className: `${feature_type.toLowerCase()}-viz-block` 
			})
		})
	})
	var dataset = new vis.DataSet(exons_dp)

	// create a data set with groups (feature types)
	var names = Object.keys(features)
	var groups = new vis.DataSet()
	for (var g = 0; g < names.length; g++) {
		groups.add({id: names[g], content: names[g]})
	}

	// DOM element where the Timeline will be attached
	var container = document.getElementById('visualization')

	// Substract dummy timestamp from start/end values so they are showed properly
	console.log(minLocation, maxLocation)
	var options = {
		format: {
			minorLabels: function(date, scale, step) {
				return date - dummyTs
			},
			majorLabels: function(date, scale, step) {
				return ''
			}
		},
		start: dummyTs,
		min: minLocation + dummyTs,
		max: maxLocation + dummyTs,
		horizontalScroll: true

	};

	// Create a Timeline
	var timeline = new vis.Timeline(container, dataset, groups, options)

	// When user clicks on the timeline, we capture the location and go to that id
	$(container).click(function (event) {
		var props = timeline.getEventProperties(event)
		var loc = props.time - dummyTs
		//scrollToAnchor('base-' + loc);
		location.hash = "#base-" + loc;
	})

	// Go to the cDNA location typed in the input box
	var goToCdna;
	$("input[name='cdna_location']").keyup(function(e) {
		loc = $(this).val()
		clearTimeout(goToCdna)
		goToCdna = setTimeout(function() {
			//scrollToAnchor('base-' + loc);
			location.hash = '#base-' + loc
			$("input[name='cdna_location']").focus()
		}, 100)
	})


	$("select#goto-region").change(function(e) {
		option = $("option:selected", this)
		loc = parseInt(option.data('start')) + 1
		//scrollToAnchor('base-' + loc);
		location.hash = '#base-' + loc
	})
}

function formatGeneSequence(sequence, features, variations, sequenceLength, basesPerRow) {
	// TODO: Refactor Sequence display as Handlebars template if possible

	// Get variations per row
	var variationsPerRow = {}
	variations.forEach(function(variation, _) {
		startRow = parseInt((variation.start - 1) / basesPerRow)
		endRow = parseInt((variation.end - 1) / basesPerRow)
		rowCount = startRow - endRow + 1
		annotatedBases = 0
		for (var i = startRow; i <= endRow; i++) {
			start = Math.max(1 + i * basesPerRow, variation.start)
			end = Math.min((i+1) * basesPerRow, variation.end)
			if (variation.operation == 'del') {
				seq = '&nbsp;'.repeat(end - start + 1)
			} else if (variation.operation == 'dup') {
				triangle = $('<span>')
				triangle.addClass('arrow-up')
				seq = $('<span>')
				seq.append(triangle)
			} else {
				seq = variation.alt.slice(annotatedBases, end - start + 1)				
			}
			annotatedBases = end - start + 1
			variationsPerRow[i] = variationsPerRow[i] || {}
			variationsPerRow[i][variation.source] = variationsPerRow[i][variation.source] || []
			variationsPerRow[i][variation.source].push({
				start: start,
				end: end,
				sequence: seq,
				variation: variation
			})
		}
	})
	console.log(variationsPerRow)

	// Sort variations per row
	Object.keys(variationsPerRow).forEach(function(row_i, _) {
		Object.keys(variationsPerRow[row_i]).forEach(function(source, _) {
			variationsPerRow[row_i][source].sort(function(a1, a2) {
				return a1.start - a2.start
			})
		})
	})

	var featuresPerRow = {}
	Object.keys(features).forEach(function(feature_type, _) {
		features[feature_type].forEach(function(feature, feature_i) {
			start = feature.location[0]
			end = feature.location[1]
			startRow = parseInt((start - 1) / basesPerRow)
			endRow = parseInt((end - 1) / basesPerRow)
			console.log(feature_type, "from row", startRow, "to", endRow)
			for (var i = startRow; i <= endRow; i++) {
				featuresPerRow[i] = featuresPerRow[i] || []
				featuresPerRow[i].push(feature_type + ' ' + (feature_i + 1))
			}
		})
	})
	console.log(featuresPerRow)

	var pos = 1
	var maxLengthDigits = sequenceLength.toString().length

	$('.col-location').css('margin-left', maxLengthDigits + 'em')

	var seqHTML = $('<p>')

	var sliceSeq, sliceTranslation, sliceVariation, row, subPos;
	subPos = 0
	row_i = 0;
	sequence.forEach(function(seqSlice, _) {
		if (seqSlice.type.toLowerCase() == 'non-cds') {
			for (base of seqSlice.sequence) {
				if ((pos-1) % basesPerRow == 0) {
					if (row !== undefined) {
						row.data('features', featuresPerRow[row_i])
						row.data('row-number', row_i)
						row.append(sliceVariation)
						row.append(sliceSeq)
						seqHTML.append(row)
						row_i++
					}
					row = $('<p>')
					rowLocation = $('<span>')

					row.addClass('seq-row')
					rowLocation.addClass('row-location')

					sliceSeq = $('<div>')
					sliceVariation = $('<div>')

					sliceSeq.addClass('seq ' + seqSlice.type.toLowerCase())
					sliceVariation.addClass('variations')

					rowLocation.html(parseInt(pos / basesPerRow) * basesPerRow)
					rowLocation.css('width', maxLengthDigits.toString() + 'em')
					sliceSeq.append(rowLocation)
					sliceVariation.css('margin-left', maxLengthDigits.toString() + 'em')

					
					// If there are variations in the row, show them
					if (variationsPerRow.hasOwnProperty(row_i)) {
						Object.keys(variationsPerRow[row_i]).forEach(function(source, _) {
							variationRow = $("<div>")
							variationRow.addClass('variation-row hidden')
							lastAnnotatedPosition = row_i * basesPerRow
							variationsPerRow[row_i][source].forEach(function(variation, _) {
								// Fill with empty spaces between variations
								start = Math.max(variation.start, row_i * basesPerRow + 1)
								console.log("Variation offset", start, lastAnnotatedPosition)
								variationRow.append("&nbsp;".repeat(start - lastAnnotatedPosition - 1))
								if (variation.variation.url) {
									variationSpan = $('<a>')
									variationSpan.attr('href', variation.variation.url)
									variationSpan.attr('target', '_blank')
								} else {
									variationSpan = $('<span>')
								}
								variationSpan.addClass('variation')
								variationSpan.addClass('variation-' + variation.variation.operation)
								variationSpan.html(variation.sequence)
								// Adding variation tooltip
								if (variation.variation.comment) {
									variationSpan.data('toggle', 'tooltip')
									variationSpan.data('placement', 'top')
									variationSpan.attr('title', variation.variation.comment)
								}
								variationRow.attr('data-source', source)
								variationRow.append(variationSpan)
								lastAnnotatedPosition = variation.end
							})
							sourceSpan = $('<span>')
							sourceSpan.addClass('variation-source')
							sourceSpan.append(source)
							variationRow.append("&nbsp;".repeat((row_i + 1) * basesPerRow - lastAnnotatedPosition + 1))
							variationRow.append(sourceSpan)
							sliceVariation.append(variationRow)
						})
					}
				}

				baseSpan = $('<span>')
				baseSpan.attr('id', 'base-' + pos)
				baseSpan.data('location', pos)
				baseSpan.addClass('base')
				// Enable tooltip
				baseSpan.data('toggle', 'tooltip')
				baseSpan.data('placement', 'top')
				baseSpan.attr('title', 'cDNA location: ' + pos)
				baseSpan.html(base)
				sliceSeq.append(baseSpan)
				pos++
			}
		}
		else if (seqSlice.type.toLowerCase() == 'cds') {
			seqSlice.triplets.forEach(function(triplet, triplet_i) {
				if (sliceTranslation == undefined) {
					sliceTranslation = $('<div>')
					sliceTranslation.addClass('translation')
					sliceTranslation.css('margin-left', maxLengthDigits.toString() + 'em')
					sliceTranslation.append("&nbsp;".repeat(pos-parseInt(pos/basesPerRow)*basesPerRow-1))
				}
				var oddEven = (triplet_i % 2 == 0) ? 'even' : 'odd'

				var tripletSeq = $('<span>')
				tripletSeq.addClass('triplet')
				var tripletTranslation = $('<span>')
				tripletTranslation.addClass('triplet ' + oddEven)

				tripletBase_i = 0
				for (base of triplet.sequence) {
					if ((pos-1) % basesPerRow == 0) {
						if (row !== undefined) {
							sliceSeq.append(tripletSeq)
							sliceTranslation.append(tripletTranslation)
							row.data('features', featuresPerRow[row_i])
							row.data('row-number', row_i)
							row.append(sliceVariation)
							row.append(sliceSeq)
							row.append(sliceTranslation)
							seqHTML.append(row)
							row_i++
						}
						row = $('<p>')
						rowLocation = $('<span>')
						sliceSeq = $('<div>')
						sliceTranslation = $('<div>')
						tripletSeq = $('<span>')
						tripletTranslation = $('<span>')

						sliceVariation = $('<div>')
						sliceVariation.addClass('variations')
						sliceVariation.css('margin-left', maxLengthDigits.toString() + 'em')

						row.addClass('seq-row')
						rowLocation.addClass('row-location')
						sliceSeq.addClass('seq ' + seqSlice.type.toLowerCase())
						sliceTranslation.addClass('translation')
						tripletSeq.addClass('triplet')
						tripletTranslation.addClass('triplet ' + oddEven)

						rowLocation.html(parseInt(pos / basesPerRow) * basesPerRow)
						rowLocation.css('width', maxLengthDigits.toString() + 'em')
						sliceSeq.append(rowLocation)
						sliceTranslation.css('margin-left', maxLengthDigits.toString() + 'em')

					// If there are variations in the row, show them
					if (variationsPerRow.hasOwnProperty(row_i)) {
						Object.keys(variationsPerRow[row_i]).forEach(function(source, _) {
							variationRow = $("<div>")
							variationRow.addClass('variation-row hidden')
							lastAnnotatedPosition = row_i * basesPerRow
							variationsPerRow[row_i][source].forEach(function(variation, _) {
								// Fill with empty spaces between variations
								start = Math.max(variation.start, row_i * basesPerRow + 1)
								console.log("Variation offset", start, lastAnnotatedPosition)
								variationRow.append("&nbsp;".repeat(start - lastAnnotatedPosition - 1))
								if (variation.variation.url) {
									variationSpan = $('<a>')
									variationSpan.attr('href', variation.variation.url)
									variationSpan.attr('target', '_blank')
								} else {
									variationSpan = $('<span>')
								}
								variationSpan.addClass('variation')
								variationSpan.addClass('variation-' + variation.variation.operation)
								variationSpan.html(variation.sequence)
								// Adding variation tooltip
								if (variation.variation.comment) {
									variationSpan.data('toggle', 'tooltip')
									variationSpan.data('placement', 'top')
									variationSpan.attr('title', variation.variation.comment)
								}
								variationRow.attr('data-source', source)
								variationRow.append(variationSpan)
								lastAnnotatedPosition = variation.end
							})
							sourceSpan = $('<span>')
							sourceSpan.addClass('variation-source')
							sourceSpan.append(source)
							variationRow.append("&nbsp;".repeat((row_i + 1) * basesPerRow - lastAnnotatedPosition + 1))
							variationRow.append(sourceSpan)
							sliceVariation.append(variationRow)
						})
					}
					}
					baseSpan = $('<span>')
					baseSpan.attr('id', 'base-' + pos)
					baseSpan.data('location', pos)
					baseSpan.addClass('base')
					// Enable tooltip
					baseSpan.data('toggle', 'tooltip')
					baseSpan.data('placement', 'top')
					baseSpan.attr('title', 'cDNA location: ' + pos)
					baseSpan.html(base)
					tripletSeq.append(baseSpan)
					tripletTranslation.append(tripletBase_i == 1 ? triplet.translation : "&nbsp;")
					pos++
					tripletBase_i++
				}
				sliceSeq.append(tripletSeq)
				sliceTranslation.append(tripletTranslation)
			})
		}
	})
	$('.sequence').html(seqHTML)
}

var currentVariations = [];
function bindVariations(popover_template) {
	$('.seq .base').mouseup(function() {
		var range = window.getSelection().getRangeAt(0)
		var $startNode = $(range.startContainer.parentNode)
		var $endNode = $(range.endContainer.parentNode)

		var startLocation = $startNode.data('location')
		var endLocation = $endNode.data('location')
		
		if (startLocation == endLocation) {
			return
		}
		var bases = range.toString()

		$('.sequence').popover('dispose')
		$('.sequence').popover({
			'placement': 'right',
			'container': 'body',
			'html': true,
			// Use Handlebars for this
			'content': popover_template({
				startLocation: startLocation,
				endLocation: endLocation,
				bases: bases.replace(/[^ATGC]/g, ''),
				singleBaseSelected: false
			})
		})
		for (var i = startLocation; i <= endLocation; i++) {
			document.getElementById('base-'+i).classList.add('selected-base')
		}
	})

	$('.seq .base').click(function(e) {
		var $node = $(this)
		$node.addClass('selected-base')
		var location = $node.data('location')
		var base = $node.text()
		$('.sequence').popover('dispose')
		$('.sequence').popover({
			'placement': 'right',
			'container': 'body',
			'html': true,
			// Use Handlebars for this
			'content': popover_template({
				startLocation: location,
				endLocation: location,
				bases: base.replace(/\s/g, ''),
				singleBaseSelected: true
			})
		})
	})

	$('body').on("change", "[name='operation']", function() {
		$option = $(this)
		operation = $option.val()
		$('.variation-sub-form').addClass('hidden')
		$('.variation-sub-form[data-operation="'+ operation +'"]').removeClass('hidden')
		$('.confirm-variation').removeClass('hidden')
	})

	$('body').on('click', '.confirm-variation', function() {
		$form = $('.variation-sub-form:visible')

		variation = getFormData($form)
		variation.operation = $form.data('operation')
		currentVariations.push(variation)
		$('#variations-form input[name="variations"]').val(JSON.stringify(currentVariations))
		$('.sequence').popover('dispose')
		console.log(currentVariations)
	})
}

function bindScroll() {
	$(".sequence").on('scroll', function(){
	   	var topRowNumber = Infinity
	   	var $topRow = $('.seq-row').eq(0)
	   	var container = $(".sequence")
		$('.seq-row').each(function() {
			var inView = checkInView($(this), container)
			if (inView && $(this).data('row-number') < topRowNumber) {
				topRowNumber = $(this).data('row-number')
				$topRow = $(this)

			}
		})
		$('.current-region').html($topRow.data('features').join(', '))
	})
	$('.sequence').scroll()
}


function bindVariationsSelect() {
	select = $('#show-variations')

	select.selectpicker({
		liveSearch: true,
		actionsBox: true,
		showTick: true
	})

	select.on('show.bs.select',function () {
		select.on('changed.bs.select', function() {
			sources = $(this).val()
			$('.variation-row').addClass('hidden')
			sources.forEach(function (source,  _) {
				$('.variation-row[data-source="'+source+'"]').removeClass('hidden')
			})
		})
	})
	select.selectpicker('refresh')
}