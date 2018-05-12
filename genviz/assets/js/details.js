var plotGeneFeatures = (features) => {
	// Dummy timestamp to use non-date values in timeline vis.js chart
	var dummyTs = +new Date()

	// Build vis.js dataset
	var exons_dp = []
	var prev_end = -1
	var row = 0
	Object.keys(features).forEach((feature_type, _) => {
		features[feature_type].forEach((feature, i) => {
			start = feature.location[0]
			end = feature.location[1]
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
	var names = Object.keys(features);
	var groups = new vis.DataSet();
	for (var g = 0; g < names.length; g++) {
		groups.add({id: names[g], content: names[g]});
	}

	// DOM element where the Timeline will be attached
	var container = document.getElementById('visualization');

	// Substract dummy timestamp from start/end values so they are showed properly
	var options = {
		format: {
			minorLabels: (date, scale, step) => {
				return date - dummyTs
			},
			majorLabels: (date, scale, step) => {
				return ''
			}
		}
	};

	// Create a Timeline
	var timeline = new vis.Timeline(container, dataset, groups, options);

	container.onclick = function (event) {
	  var props = timeline.getEventProperties(event)
	  var location = props.time - dummyTs
	  console.log(location)
	}
}