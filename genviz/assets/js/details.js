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

var formatGeneSequence = function(sequence, annotations, sequenceLength, basesPerRow) {
	// TODO: Refactor URGENT!!!
	var pos = 1
	var maxLengthDigits = sequenceLength.toString().length

	$('.col-location').css('margin-left', maxLengthDigits + 'em')

	var seqHTML = $('<p>')

	var sliceSeq, sliceTranslation, sliceAnnotation, row, subPos;
	subPos = 0
	sequence.forEach(function(seqSlice, _) {
		if (seqSlice.type.toLowerCase() == 'non-cds') {
			for (base of seqSlice.sequence) {
				if ((pos-1) % basesPerRow == 0) {
					if (row !== undefined) {
						row.append(sliceSeq)
						row.append(sliceAnnotation)
						seqHTML.append(row)
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
				}
				annotationSpan = $('<span>')
				annotationSpan.addClass('annotation')
				is_annotated = false
				annotations.forEach(function(annotation, _) {
					if (annotation.operation == 'ins' && annotation.after == pos-1) {
						annotationSpan.addClass("annotation-ins")
						annotationSpan.html(annotation.sequence)
						is_annotated = true
					}
					else if (pos >= annotation.start && pos <= annotation.end) {
						if (annotation.operation == 'del') {
							annotationSpan.addClass("annotation-del")
							annotationSpan.html("&nbsp;")
							is_annotated = true
						}
						else if (annotation.operation == 'sub') {
							annotationSpan.addClass("annotation-sub")
							annotationSpan.html(annotation.sequence[subPos])
							if (++subPos == annotation.sequence.length) {
								subPos = 0
							}
							is_annotated = true
						}
					}
				})
				if (is_annotated) {
					sliceAnnotation.append(annotationSpan)
				} else {
					sliceAnnotation.append("&nbsp;")
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
function bindAnnotations() {
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
			'container': '.sequence',
			'html': true,
			// Use Handlebars for this
			'content': `<div class="form-group row">
    <label for="Operation" class="col-sm-12 col-form-label">Operation</label>
    <div class="col-sm-4">
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-del" value="del">
		  <label class="form-check-label" for="Del">Del</label>
		</div>
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-insdel" value="insdel">
		  <label class="form-check-label" for="InsDel">Insdel</label>
		</div>
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-sub" value="sub">
		  <label class="form-check-label" for="Sub">Sub</label>
		</div>
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-dup" value="dup">
		  <label class="form-check-label" for="Dup">Dup</label>
		</div>
    </div>
    <div class="col-sm-8">
       	<form class="hidden annotation-sub-form" data-operation='del'>
    		<span>Delete bases <span class='popover-bases'>(${bases.replace(/\s/g, '')})</span> from position ${startLocation} to ${endLocation}</span>
    		<input type="hidden" name="start" value="${startLocation}" />
    		<input type="hidden" name="end" value="${endLocation}" />
    	</form>
       	<form class="hidden annotation-sub-form" data-operation='insdel'>
    		<label>Delete bases <span class='popover-bases'>(${bases.replace(/\s/g, '')})</span> from position ${startLocation} to ${endLocation} and insert:</label>
    		<input type="text" class="form-control" name="sequence" />
    		<input type="hidden" name="start" value="${startLocation}" />
    		<input type="hidden" name="end" value="${endLocation}" />
    	</form>
       	<form class="hidden annotation-sub-form" data-operation='sub'>
    		<label>Substitute bases <span class='popover-bases'>(${bases.replace(/\s/g, '')})</span> from position ${startLocation} to ${endLocation} with:</label>
    		<input type="text" class="form-control" maxlength=${endLocation - startLocation + 1} name="sequence" />
    		<input type="hidden" name="start" value="${startLocation}" />
    		<input type="hidden" name="end" value="${endLocation}" />
    	</form>
       	<form class="hidden annotation-sub-form" data-operation='dup'>
    		<span>Duplicate bases <span class='popover-bases'>(${bases.replace(/\s/g, '')})</span> from position ${startLocation} to ${endLocation}</span>
    		<input type="hidden" name="start" value="${startLocation}" />
    		<input type="hidden" name="end" value="${endLocation}" />
    	</form>
    	<button type="button" class="btn btn-primary confirm-annotation hidden">Confirm</button>
    </div>
</div>
`
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
			'container': '.sequence',
			'html': true,
			// Use Handlebars for this
			'content': `<div class="form-group row">
    <label for="Operation" class="col-sm-12 col-form-label">Operation</label>
    <div class="col-sm-4">
      	<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-ins" value="ins">
		  <label class="form-check-label" for="Ins">Ins</label>
		</div>
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-del" value="del">
		  <label class="form-check-label" for="Del">Del</label>
		</div>
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-insdel" value="insdel">
		  <label class="form-check-label" for="InsDel">InsDel</label>
		</div>
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-sub" value="sub">
		  <label class="form-check-label" for="Sub">Sub</label>
		</div>
		<div class="form-check">
		  <input class="form-check-input" type="radio" name="operation" id="radio-dup" value="dup">
		  <label class="form-check-label" for="Dup">Dup</label>
		</div>
    </div>
    <div class="col-sm-8">
    	<form class="hidden annotation-sub-form" data-operation='ins'>
    		<label>Insert after position ${location}:</label>
    		<input type="text" class="form-control" name="sequence" />
    		<input type="hidden" name="after" value="${location}" />
    	</form>
       	<form class="hidden annotation-sub-form" data-operation='del'>
    		<span>Delete base (${base}) at position ${location}</span>
    		<input type="hidden" name="start" value="${location}" />
    		<input type="hidden" name="end" value="${location}" />
    	</form>
       	<form class="hidden annotation-sub-form" data-operation='insdel'>
    		<label>Delete base (${base}) at position ${location} and insert:</label>
    		<input type="text" class="form-control" name="sequence" />
    		<input type="hidden" name="start" value="${location}" />
    		<input type="hidden" name="end" value="${location}" />
    	</form>
       	<form class="hidden annotation-sub-form" data-operation='sub'>
    		<label>Substitute base (${base}) at position ${location} with:</label>
    		<input type="text" class="form-control" maxlength=1 name="sequence" />
    		<input type="hidden" name="start" value="${location}" />
    		<input type="hidden" name="end" value="${location}" />
    	</form>
       	<form class="hidden annotation-sub-form" data-operation='dup'>
    		<span>Duplicate base (${base}) at position ${location}</span>
    		<input type="hidden" name="start" value="${location}" />
    		<input type="hidden" name="end" value="${location}" />
    	</form>
    	<button type="button" class="btn btn-primary confirm-annotation hidden">Confirm</button>
    </form>
</div>
`
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
