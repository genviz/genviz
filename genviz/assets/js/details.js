function scrollToAnchor(id){
    var tag = $("#"+id);
    $('html,body').animate({scrollTop: tag.offset().top});
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
		scrollToAnchor('base-' + loc);
		location.hash = "#base-" + loc;
	})

	// Go to the cDNA location typed in the input box
	var goToCdna;
	$("input[name='cdna_location']").keyup(function(e) {
		loc = $(this).val()
		clearTimeout(goToCdna)
		goToCdna = setTimeout(function() {
			scrollToAnchor('base-' + loc);
			location.hash = '#base-' + loc
			$("input[name='cdna_location']").focus()
		}, 100)
	})
}

function formatGeneSequence(sequence, annotations, sequenceLength, basesPerRow) {
	// TODO: Refactor URGENT!!!

	// Get annotations per row
	var annotationsPerRow = {}
	annotations.forEach(function(annotation, _) {
		startRow = parseInt((annotation.start - 1) / basesPerRow)
		endRow = parseInt((annotation.end - 1) / basesPerRow)
		rowCount = startRow - endRow + 1
		annotatedBases = 0
		for (var i = startRow; i <= endRow; i++) {
			start = Math.max(1 + i * basesPerRow, annotation.start)
			end = Math.min((i+1) * basesPerRow, annotation.end)
			if (annotation.operation == 'del') {
				seq = '&nbsp;'.repeat(end - start + 1)
			} else {
				seq = annotation.sequence.slice(annotatedBases, end - start + 1)				
			}
			annotatedBases = end - start + 1
			annotationsPerRow[i] = annotationsPerRow[i] || {}
			annotationsPerRow[i][annotation.source] = annotationsPerRow[i][annotation.source] || []
			annotationsPerRow[i][annotation.source].push({
				start: start,
				end: end,
				sequence: seq,
				annotation: annotation
			})
		}
	})
	// Sort annotations per row
	Object.keys(annotationsPerRow).forEach(function(row_i, _) {
		Object.keys(annotationsPerRow[row_i]).forEach(function(source, _) {
			annotationsPerRow[row_i][source].sort(function(a1, a2) {
				return a1.start - a2.start
			})
		})
	})
	console.log(annotationsPerRow)

	var pos = 1
	var maxLengthDigits = sequenceLength.toString().length

	$('.col-location').css('margin-left', maxLengthDigits + 'em')

	var seqHTML = $('<p>')

	var sliceSeq, sliceTranslation, sliceAnnotation, row, subPos;
	subPos = 0
	row_i = 0;
	sequence.forEach(function(seqSlice, _) {
		if (seqSlice.type.toLowerCase() == 'non-cds') {
			for (base of seqSlice.sequence) {
				if ((pos-1) % basesPerRow == 0) {
					if (row !== undefined) {
						row.append(sliceAnnotation)
						row.append(sliceSeq)
						seqHTML.append(row)
						row_i++
					}
					row = $('<p>')
					rowLocation = $('<span>')

					row.addClass('seq-row')
					rowLocation.addClass('row-location')

					sliceSeq = $('<div>')
					sliceAnnotation = $('<div>')

					sliceSeq.addClass('seq ' + seqSlice.type.toLowerCase())
					sliceAnnotation.addClass('annotations')

					rowLocation.html(parseInt(pos / basesPerRow) * basesPerRow)
					rowLocation.css('width', maxLengthDigits.toString() + 'em')
					sliceSeq.append(rowLocation)
					sliceAnnotation.css('margin-left', maxLengthDigits.toString() + 'em')

					
					// If there are annotations in the row, show them
					if (annotationsPerRow.hasOwnProperty(row_i)) {
						Object.keys(annotationsPerRow[row_i]).forEach(function(source, _) {
							annotationRow = $("<div class='annotation-row'>")
							lastAnnotatedPosition = row_i * basesPerRow
							annotationsPerRow[row_i][source].forEach(function(annotation, _) {
								// Fill with empty spaces between annotations
								annotationRow.append("&nbsp;".repeat(annotation.start - lastAnnotatedPosition - 1))
								annotationSpan = $('<span>')
								annotationSpan.addClass('annotation')
								annotationSpan.addClass('annotation-' + annotation.annotation.operation)
								annotationSpan.html(annotation.sequence)
								// Adding annotation tooltip
								if (annotation.annotation.comment) {
									annotationSpan.data('toggle', 'tooltip')
									annotationSpan.data('placement', 'top')
									annotationSpan.attr('title', 'Comment: ' + annotation.annotation.comment)
								}
								annotationRow.append(annotationSpan)
								lastAnnotatedPosition = annotation.end
							})
							sourceSpan = $('<span>')
							sourceSpan.addClass('annotation-source')
							sourceSpan.append(source)
							annotationRow.append("&nbsp;".repeat((row_i + 1) * basesPerRow - lastAnnotatedPosition + 1))
							annotationRow.append(sourceSpan)
							sliceAnnotation.append(annotationRow)
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
							row.append(sliceSeq)
							row.append(sliceTranslation)
							seqHTML.append(row)
						}
						row = $('<p>')
						rowLocation = $('<span>')
						sliceSeq = $('<div>')
						sliceTranslation = $('<div>')
						tripletSeq = $('<span>')
						tripletTranslation = $('<span>')

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

var currentAnnotations = [];
function bindAnnotations(popover_template) {
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
		$('.annotation-sub-form').addClass('hidden')
		$('.annotation-sub-form[data-operation="'+ operation +'"]').removeClass('hidden')
		$('.confirm-annotation').removeClass('hidden')
	})

	$('body').on('click', '.confirm-annotation', function() {
		$form = $('.annotation-sub-form:visible')

		annotation = getFormData($form)
		annotation.operation = $form.data('operation')
		currentAnnotations.push(annotation)
		$('#annotations-form input[name="annotations"]').val(JSON.stringify(currentAnnotations))
		$('.sequence').popover('dispose')
		console.log(currentAnnotations)
	})
}
