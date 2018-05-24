function scrollToAnchor(id){
    var tag = $("#"+id);
    $('html,body').animate({scrollTop: tag.offset().top});
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

	container.onclick = function (event) {
	  var props = timeline.getEventProperties(event)
	  var loc = props.time - dummyTs
	  scrollToAnchor('base-' + loc);
	  location.hash = "#base-" + loc;
	}

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

var formatGeneSequence = function(sequence, sequenceLength, basesPerRow) {
	// TODO: Refactor
	var pos = 0
	var maxLengthDigits = sequenceLength.toString().length

	$('.col-location').css('margin-left', maxLengthDigits + 'em')

	var seqHTML = $('<p>')

	var sliceSeq, sliceTranslation, row;
	sequence.forEach(function(seqSlice, _) {
		if (seqSlice.type.toLowerCase() == 'non-cds') {
			for (base of seqSlice.sequence) {
				if (pos % basesPerRow == 0) {
					if (row !== undefined) {
						row.append(sliceSeq)
						row.append(sliceTranslation)
						seqHTML.append(row)
					}
					row = $('<p>')
					rowLocation = $('<span>')

					row.addClass('seq-row')
					rowLocation.addClass('row-location')

					sliceSeq = $('<div>')
					sliceTranslation = $('<div>')

					sliceSeq.addClass('seq ' + seqSlice.type.toLowerCase())
					sliceTranslation.addClass('translation')

					rowLocation.html(parseInt(pos / basesPerRow) * basesPerRow)
					rowLocation.css('width', maxLengthDigits.toString() + 'em')
					sliceSeq.append(rowLocation)
					sliceTranslation.css('margin-left', maxLengthDigits.toString() + 'em')
				}
				baseSpan = $('<span>')
				baseSpan.attr('id', 'base-' + pos)
				baseSpan.html(base)
				sliceSeq.append(baseSpan)
				sliceTranslation.append("&nbsp;")
				pos++
			}
		}
		else if (seqSlice.type.toLowerCase() == 'cds') {
			seqSlice.triplets.forEach(function(triplet, triplet_i) {
				var oddEven = (triplet_i % 2 == 0) ? 'even' : 'odd'

				var tripletSeq = $('<span>')
				tripletSeq.addClass('triplet')
				var tripletTranslation = $('<span>')
				tripletTranslation.addClass('triplet ' + oddEven)

				tripletBase_i = 0
				for (base of triplet.sequence) {
					if (pos % basesPerRow == 0) {
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