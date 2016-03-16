function stickyScatter() 
{
    $('#scatterContainer').addClass('stick');
	var top = Math.max($(window).scrollTop() + 10, $('#confidences').offset().top);
    var left = $('#confidences').offset().left + $('#confidences').width() + 25;
	$('#scatterContainer').offset({ left: left, top: top });
}

function arrayUnique(arr) 
{
    return arr.reduce(function(p, c) 
    {
        if (p.indexOf(c) < 0) p.push(c);
        return p;
    }, []);
}

function setActiveExclusively(elem)
{
	$('.active').removeClass('active');
	elem.addClass('active');
}

function shadeColor(color, percent) 
{   
    var f=parseInt(color.slice(1),16),t=percent<0?0:255,p=percent<0?percent*-1:percent,R=f>>16,G=f>>8&0x00FF,B=f&0x0000FF;
    return "#"+(0x1000000+(Math.round((t-R)*p)+R)*0x10000+(Math.round((t-G)*p)+G)*0x100+(Math.round((t-B)*p)+B)).toString(16).slice(1);
}

function blendColors(c0, c1, p) 
{
    var f=parseInt(c0.slice(1),16),t=parseInt(c1.slice(1),16),R1=f>>16,G1=f>>8&0x00FF,B1=f&0x0000FF,R2=t>>16,G2=t>>8&0x00FF,B2=t&0x0000FF;
    return "#"+(0x1000000+(Math.round((R2-R1)*p)+R1)*0x10000+(Math.round((G2-G1)*p)+G1)*0x100+(Math.round((B2-B1)*p)+B1)).toString(16).slice(1);
}

function onlyUnique(value, index, self) 
{ 
    return self.indexOf(value) === index;
}

function removeIndexes(arr, rmIdx)
{
	rmIdx.sort(function(a,b){return a - b});
	rmIdx.reverse();
	var narr = arr.slice();

	for (var i = 0; i < rmIdx.length; i++)
	{
		narr.splice(rmIdx[i], 1);
	}

	return narr;
}

function cellBarChart(containerCell, confidences, maxK)
{
	// convert confidences into data array
	var data = new Array(maxK);
	for (var i = 1; i <= maxK; i++) 
	{
		data[i-1] = {};
		data[i-1].numClusters = i;
		if (i in confidences)
		{
			data[i-1].confidence = confidences[i];
		} else
		{
			data[i-1].confidence = 0;
		}
	}

	var margin = {top: 10, right: 20, bottom: 25, left: 20},
	width = 160 - margin.left - margin.right,
	height = 60 - margin.top - margin.bottom;

	var x = d3.scale.ordinal().rangeRoundBands([0, width], .1);
	var y = d3.scale.linear().range([height, 0]);

	var xAxis = d3.svg.axis().outerTickSize(0).scale(x).orient("bottom");
	var yAxis = d3.svg.axis().outerTickSize(0).scale(y).orient("left").ticks(0);

	var svg = d3.select($(containerCell).get(0)).append("svg")
		.attr("width", width + margin.left + margin.right)
		.attr("height", height + margin.top + margin.bottom)
		.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	x.domain(data.map(function(d) { return d.numClusters; }));
	y.domain([0, 1]);

	svg.append("rect")
	   .attr("x", 0)
	   .attr("y", 0)
	   .attr("width", width)
	   .attr("height", height)
	   .style("fill", "none")
  	   .style("stroke", "#aaa")
       .style("stroke-width", "1");	  

  	svg.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(0," + height + ")")
		.call(xAxis)
		.append("text")
		.attr("y", 25)
		.text("# clusters");

  	svg.append("g")
		.attr("class", "y axis")
		.call(yAxis)
		.append("text")
		.attr("transform", "rotate(-90)")
		.attr("y", -5)
		.attr("x", 5)
		.style("text-anchor", "end")
		.text("confidence");	

    tip = d3.tip().attr('class', 'd3-tip').offset([-10, 0]).html(function(d, i) { return 'p(' + data[i].numClusters + ') = ' + data[i].confidence; });
    svg.call(tip);

	svg.selectAll(".bar")
		.data(data)
		.enter().append("rect")
		.attr("class", "bar")
		.attr("x", function(d) { return x(d.numClusters) + 2; })
		.attr("width", x.rangeBand()-4)
		.attr("y", function(d) { return y(d.confidence); })
		.attr("height", function(d) { return height - y(d.confidence); })
		.on('mouseover', tip.show)
		.on('mouseout', tip.hide);

}

function buildConfidenceTable(results)
{
	if (typeof results == 'undefined' || Object.keys(results) == 0)
	{
		return;
	}

	var krakenEnabled = "krakenLabels" in results[Object.keys(results)[0]];

	$('#confidences').append('<tr><th>ID</th><th>Contamination<br/>status</th><th>Sample</th><th>CC</th><th>Dip</th>' + (krakenEnabled ? '<th>Kraken</th>' : '') + '</tr>');

	for (var i in results)
	{
		var kraken = '';
		if (krakenEnabled)
		{
			var krakenNumUnknown = 0;
			$.each(results[i].krakenLabels, function (idx, val) 
			{
				if (val === "unknown")
				{
					krakenNumUnknown++;
				}
			});
			var krakenUnique = results[i].krakenLabels.filter(onlyUnique);
			var numSpecies = krakenUnique.length;
			kraken = krakenNumUnknown > 0 ? '&ge; ' + numSpecies : '' + numSpecies;
		}		

		var contProbCC = 0;
		var contProbDip = 0;
		for (var j = 0; j < results[i].bootstraps.length; j++)
		{
			contProbCC += results[i].bootstraps[j].hasSeparatedComponents ? 1 : 0;
			contProbDip += results[i].bootstraps[j].isMultiModal ? 1 : 0;
		}
		contProbCC /= results[i].bootstraps.length;
		contProbDip /= results[i].bootstraps.length;

		var status = 'warning';
		if (contProbCC < 0.25 && contProbDip < 0.25)
		{
			status = 'clean';
		} else if (contProbCC > 0.75 || contProbDip > 0.75)
		{
			status = 'contaminated';
		}	

		$('#confidences').append('<tr>' +
			'<td>' + results[i].id + '</td>' +
			'<td class="' + status + '">&nbsp;</td>' + 
			'<td class="selectable dataConf">' + i + '</td>' +
			'<td class="selectable ccConf">p = ' + contProbCC.toFixed(2) + '</td>' +
			'<td class="selectable dipConf">p = ' + contProbDip.toFixed(2) + '</td>' +
			(krakenEnabled ? '<td class="selectable kraken"><span class="number">' + kraken + '</span><br/>species</td>' : '') +
			'</tr>');

	}	
}

function dimMax(data, dim)
{
	return data.reduce(function(max, arr) 
	{
	    return max >= arr[dim] ? max : arr[dim];
	}, -Infinity);
}

function dimMin(data, dim)
{
	return data.reduce(function(min, arr) 
	{
	    return min < arr[dim] ? min : arr[dim];
	}, Infinity);
}

function showData(dataMat, labels, tooltips, greyedOut, width, height, padding)
{	
	var svg = d3.select("#scatter").html("").append("svg").attr("width", width).attr("height", height);

	var doGrey = typeof greyedOut != 'undefined' && greyedOut.length == labels.length;

	svg.append("rect")
	   .attr("x", 0)
	   .attr("y", 0)
	   .attr("width", width)
	   .attr("height", height)
	   .style("fill", "none")
  	   .style("stroke", "lightgrey")
       .style("stroke-width", "1");	   

	minX = dimMin(dataMat, 0);
	maxX = dimMax(dataMat, 0);
	minY = dimMin(dataMat, 1);
	maxY = dimMax(dataMat, 1);

	var xScale = d3.scale.linear()
			.domain([minX, maxX])
			.range([padding/2, width-padding/2]);
	var yScale = d3.scale.linear()
             .domain([minY, maxY])
             .range([padding/2, height-padding/2]);

    tip = d3.tip().attr('class', 'd3-tip').offset([-10, 0]).html(function(d, i) { return tooltips[i]; });
    svg.call(tip);

	svg.selectAll("circle")
			.data(dataMat)
			.enter()
			.append("circle")
			.attr("cx", function(d) { return xScale(d[0]); })
       		.attr("cy", function(d) { return yScale(d[1]); })
       		.attr("r", function(d) {return 3; })
       		.style("fill", function(d, i) {return (doGrey && greyedOut[i]) ? '#b3b3b3' : colors[labels[i] % colors.length]; })
			.on('mouseover', tip.show)
			.on('mouseout', tip.hide);	
}

function showVisualization()
{
	var x = results[selectedFasta];

	var clustAnaResult;
	if (selectedData === 'oneshot')
	{
		clustAnaResult = x.oneshot;
	} else
	{
		clustAnaResult = x.bootstraps[parseInt(selectedData)];
	}

	var dataMat = clustAnaResult[selectedReduction];

	var labels;
	var greyedOut = new Array();
	var outliers = new Array();
	if (selectedLabels === 'fasta')
	{
		// labels = numericLabels(x.fastaLabels);
		labels = Array.apply(null, Array(x.fastaLabels.length)).map(Number.prototype.valueOf, -1); // no labels / black color
	} else if(selectedLabels === 'cc')
	{
		labels = clustAnaResult.clustCC.labels;	
		outliers = clustAnaResult.clustCC.outlierClusters;	
	} else if(selectedLabels === 'dip')
	{
		if (selectedReduction === 'dataPca')
		{
			if (!selectedNumClusters)
			{
				selectedNumClusters = clustAnaResult.numClustPca;
				$('#numClusters' + selectedNumClusters).prop('checked', true);
			}
			labels = clustAnaResult.clustsPca[selectedNumClusters-1].labels;
			outliers = clustAnaResult.clustsPca[selectedNumClusters-1].outlierClusters;	
		} else if (selectedReduction === 'dataSne')
		{
			if (!selectedNumClusters)
			{
				selectedNumClusters = clustAnaResult.numClustSne;
				$('#numClusters' + selectedNumClusters).prop('checked', true);
			}			
			labels = clustAnaResult.clustsSne[selectedNumClusters-1].labels;
			outliers = clustAnaResult.clustsSne[selectedNumClusters-1].outlierClusters;	
		}
	} else if(selectedLabels === 'kraken')
	{
		labels = numericLabels(x.krakenLabels);
		for (var i = 0; i < x.krakenLabels.length; i++)
		{
			if (x.krakenLabels[i] === 'unknown')
			{
				greyedOut.push(true);
			} else
			{
				greyedOut.push(false);
			}
		}
	}	

	var tooltips = results[selectedFasta].fastaLabels;
	if(selectedLabels === 'kraken')
	{
		tooltips = x.krakenLabels;
	} else
	{
		if (selectedData !== 'oneshot') // align bootstrap tooltips
		{
			tooltips = bootstrapLabels(tooltips, clustAnaResult.bootstrapIndexes);
		}
	}


	if (!showOutliers && typeof outliers != 'undefined' && outliers.length > 0)
	{		
		var n = dataMat.length;
		var rmIdx = new Array();

		for (var i = 0; i < n; i++)
		{
			if ($.inArray(labels[i], outliers) >= 0)
			{
				rmIdx.push(i);
			}			
		}
		rmIdx = arrayUnique(rmIdx);

		dataMat = removeIndexes(dataMat, rmIdx);
		labels = removeIndexes(labels, rmIdx);
		tooltips = removeIndexes(tooltips, rmIdx);
	}

	showData(dataMat, labels, tooltips, greyedOut, width, height, padding);
	updateExport(labels);
}

function updateBootStrapsSelect()
{
	$('#bootstraps').val('oneshot');
	selectedData = 'oneshot';

	if (selectedLabels === 'kraken')
	{
		$('#bootstraps').prop("disabled", true);
		return;
	}

	$('#bootstraps').prop("disabled", false);

	var n = results[selectedFasta].bootstraps.length;

	$('#bootstraps').html('');
	$('#bootstraps').append('<option value="oneshot">One shot</option>');
	for (var i = 0; i < n; i++) 
	{
		$('#bootstraps').append($('<option>', {
		    value: i,
		    text: 'Bootstrap ' + (i+1)
		}));		
	}

}

function updateExport(labels)
{
	if (selectedLabels === 'fasta' || selectedData !== 'oneshot')
	{
		$('#export').hide();
		return;
	}

	var uniqueLabels = arrayUnique(labels).sort();

	$('#exportColors').html('');
	for (var i in uniqueLabels)
	{
		var lbl = uniqueLabels[i];
		var color = colors[lbl % colors.length];
		if (selectedLabels !== 'dip')
		{
			var href = 'export/' + results[selectedFasta].id + "-" + selectedLabels + "-" + lbl + ".fasta";
		} else
		{
			var href = 'export/' + results[selectedFasta].id + "-dip-" + selectedReduction + "-" + selectedNumClusters + "-" + lbl + ".fasta";
		}
		$('#exportColors').append('<span><a style="color: ' + color + ';" href="' + href + '">&#x25CF;</a></span>&nbsp;');
	}

	$('#export').show();
}

function numericLabels(labels)
{
	var mp = new Array();
	k = 0;
	for (var i in labels) 
	{
		var key = labels[i];
		if (!(key in mp))
		{
			mp[key] = k;
			k++;
		}
	}

	var result = Array(labels.length);
	for (var i in labels) 
	{
		result[i] = mp[labels[i]];
	}
	return result;
}

function bootstrapLabels(labelsOneshot, bootstrapIndexes)
{
	result = Array(bootstrapIndexes.length);
	for (var i in bootstrapIndexes)
	{
		result[i] = labelsOneshot[bootstrapIndexes[i]];
	}	
	return result;
}
