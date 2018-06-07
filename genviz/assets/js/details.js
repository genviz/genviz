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

    return isTotal
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

function splitVariationsInRows(variations, basesPerRow) {
	var variationsPerRow = {}
	variations.forEach(function(variation, _) {
		startRow = parseInt((variation.start - 1) / basesPerRow)
		endRow = parseInt((variation.end - 1) / basesPerRow)
		rowCount = startRow - endRow + 1
		annotatedBases = 0
		for (var i = startRow; i <= endRow; i++) {
			start = Math.max(1 + i * basesPerRow, variation.start)
			end = Math.min((i+1) * basesPerRow, variation.end)
			if (variation.operation != 'del' && variation.operation != 'dup') {
				seq = variation.alt.substr(annotatedBases, end - start + 1)				
			} else {
				seq = variation.alt ? variation.alt.slice(annotatedBases, end - start + 1) : ""				
			}
			annotatedBases = end - start + 1
			variationsPerRow[i] = variationsPerRow[i] || {}
			variationsPerRow[i][variation.source] = variationsPerRow[i][variation.source] || []
			variationsPerRow[i][variation.source].push({
				start: start,
				end: end,
				sequence: seq,
				operation: variation.operation,
				changeLength: end - start + 1,
				variation: variation
			})
		}
	})
	console.log(variationsPerRow)

	Object.keys(variationsPerRow).forEach(function(row_i, _) {
		Object.keys(variationsPerRow[row_i]).forEach(function(source, _) {
			// Sort variations per row
			variationsPerRow[row_i][source].sort(function(a1, a2) {
				return a1.start - a2.start
			})
			// Compute spaces between variations and offset from the beginning of the row
			lastAnnotatedPosition = row_i * basesPerRow
			variationsPerRow[row_i][source].forEach(function(variation, i, vars) {
				vars[i].offset = variation.start - lastAnnotatedPosition - 1
				vars[i].changeLength = variation.end - variation.start + 1
				lastAnnotatedPosition = variation.end
			})
		})
	})
	return variationsPerRow
}

function getFeaturesPerRow(features, basesPerRow) {
	// Get features in each row
	var featuresPerRow = {}
	var translationPerRow = {}
	Object.keys(features).forEach(function(feature_type, _) {
		features[feature_type].forEach(function(feature, feature_i) {
			startFeature = feature.location[0]
			endFeature = feature.location[1]
			startRow = parseInt((startFeature - 1) / basesPerRow)
			endRow = parseInt((endFeature - 1) / basesPerRow)
			translationInserted = 0
			for (var i = startRow; i <= endRow; i++) {
				featuresPerRow[i] = featuresPerRow[i] || {}
				featuresPerRow[i].features = featuresPerRow[i].features || []
				featuresPerRow[i].features.push(feature)
				if (feature_type.toLowerCase() == 'cds') {
					start = Math.max(1 + i * basesPerRow, startFeature)
					end = Math.min((i+1) * basesPerRow, endFeature)
					offset = start - i * basesPerRow - 1
					translation = features[feature_type][feature_i].qualifiers.translation[0]
					translationToInsert = parseInt((end - start + 1) / 3)
					translationInRow = translation.substr(translationInserted, translationToInsert)
					translationInserted += translationToInsert
					featuresPerRow[i].translation = {
						translation: translationInRow,
						offset: offset
					}
				}
			}
		})
	})
	return featuresPerRow
}

function formatGeneSequence(sequence, features, variations, sequenceLength, basesPerRow, template) {
	// Get variations per row
	var variationsPerRow = splitVariationsInRows(variations, basesPerRow)
	var featuresPerRow = getFeaturesPerRow(features, basesPerRow)
	
	var maxLengthDigits = sequenceLength.toString().length

	var rowCount = parseInt(sequenceLength / basesPerRow)

	var rows = {}
	for (var row_i = 0; row_i <= rowCount; row_i++) {
		var seq = sequence.substr(row_i * basesPerRow, basesPerRow)
		var rowCDS = featuresPerRow[row_i].features.findIndex(function(f) { return f.type.toLowerCase() == 'cds' })

		rows[row_i] = {
			sequence: seq.split('').map(function(base, i) {
				pos = basesPerRow * row_i + i + 1
				return {
					position: pos,
					base: base,
					isCDS: rowCDS !== -1 && pos >= featuresPerRow[row_i].features[rowCDS][0] && pos <= featuresPerRow[row_i].features[rowCDS][1]
				}
			}),
			variations: variationsPerRow[row_i],
			translation: featuresPerRow[row_i].translation
		}
	}
	console.log(rows)

	// Render Handlebars template
	$('.sequence-container').html(template({
		rows: rows,
		maxLengthDigits,
		rowCount: rowCount,
		basesPerRow: basesPerRow
	}))
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