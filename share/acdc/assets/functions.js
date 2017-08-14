var isFileSaverSupported = false;
try 
{
    var isFileSaverSupported = !!new Blob;
} catch (e) 
{
}
if(!isFileSaverSupported)
{
	alert('HTML5 Blobs are not supported and fasta export will not work. Please upgrade your browser.');
}

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

function calculateStar(centerX, centerY, arms, outerRadius, innerRadius)
{
   var results = "";

   var angle = Math.PI / arms;

   for (var i = 0; i < 2 * arms; i++)
   {
      // Use outer or inner radius depending on what iteration we are in.
      var r = (i & 1) == 0 ? outerRadius : innerRadius;

      var currX = centerX + Math.cos(i * angle) * r;
      var currY = centerY + Math.sin(i * angle) * r;

      // Our first time we simply append the coordinates, subsequet times
      // we append a ", " to distinguish each coordinate pair.
      if (i == 0)
      {
         results = currX + "," + currY;
      }
      else
      {
         results += ", " + currX + "," + currY;
      }
   }

   return results;
}

function buildConfidenceTable()
{
	if (typeof results == 'undefined' || Object.keys(results) == 0)
	{
		return;
	}

	var krakenEnabled = results[Object.keys(results)[0]].kraken.enabled;

	$('#confidences').append('<tr><th>ID</th><th>Contamination<br/>state</th><th>Sample</th><th>CC</th><th>Dip</th>' + (krakenEnabled ? '<th>Kraken</th>' : '') + '</tr>');

	for (var i in results)
	{
		var kraken = '';
		if (krakenEnabled)
		{
			var krakenNumUnknown = 0;
			$.each(results[i].kraken.classification, function (idx, val)
			{
				if (val === "unknown")
				{
					krakenNumUnknown++;
				}
			});
			var krakenUnique = results[i].kraken.classification.filter(onlyUnique);
			var numSpecies = krakenUnique.length;
			kraken = krakenNumUnknown > 0 ? '&ge; ' + numSpecies : '' + numSpecies;
		}

		$('#confidences').append('<tr>' +
			'<td>' + results[i].id + '</td>' +
			'<td class="' + results[i].contamination_state + '">&nbsp;</td>' +
			'<td class="selectable dataConf">' + i + '</td>' +
			'<td class="selectable ccConf">conf = ' + results[i].confidence_cc.toFixed(2) + '</td>' +
			'<td class="selectable dipConf">conf = ' + results[i].confidence_dip.toFixed(2) + '</td>' +
			(krakenEnabled ? '<td class="selectable kraken"><span class="number">' + kraken + '</span> species</td>' : '') +
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

function showData(id, dataMat, labels, tooltips, greyedOut, sixteenS, width, height, padding)
{
	var svg = d3.select("#scatter").html("").append("svg").attr("width", width).attr("height", height);

	var doGrey = typeof greyedOut != 'undefined' && greyedOut.length == labels.length;
	var do16S = typeof sixteenS != 'undefined' && sixteenS.length == labels.length;

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

	var shapes = svg.selectAll("circle").data(dataMat).enter();

	shapes.append("circle")
		  .attr("cx", function(d) { return xScale(d[0]); })
		  .attr("cy", function(d) { return yScale(d[1]); })
		  .attr("r", function(d) {return 3; })
		  .style("fill", function(d, i) {return (doGrey && greyedOut[i]) ? greyColor : colors[labels[i] % colors.length]; })
		  .on('mouseover', tip.show)
		  .on('mouseout', tip.hide);

	if (do16S)
	{
		shapes.append("svg:polygon")
			  .select(function(d, i){ return sixteenS[i] ? this : null; })
			  .attr("id", "star_1")
			  .attr("visibility", "visible")
			  .attr("points", function(d) { return calculateStar(xScale(d[0]), yScale(d[1]), 8, 50, 6); })
			  .style("stroke", highlightColor)
			  .style("stroke-width", 1)
			  .style("fill", highlightColor)
			  .style("fill-opacity", .5)
			  .on('mouseover', tip.show)
			  .on('mouseout', tip.hide)
			  .on('click', function(d, i) { export16S(i); });
	}

}

function labelsPerPoint(x, labelsPerContig)
{
	var map = new Array();
	for (var i = 0; i < x.fasta_stats.included_contigs.length; i++) 
	{
		c = x.fasta_stats.included_contigs[i];
		map[c] = labelsPerContig[i];
	}	

	var res = new Array();
	for (var i = 0; i < x.visualizations.contig_labels.length; i++) 
	{
		contigIdx = x.visualizations.contig_labels[i]
		c = x.fasta_stats.included_contigs[contigIdx]
		res[i] = map[c]
	}

	return res;
}

function showVisualization()
{
	var x = results[selectedFasta];

	var dataMat_;
	if(selectedReduction === 'dataSne')
	{
		dataMat_ = x.visualizations.sne;
	}
	if(selectedReduction === 'dataPca')
	{
		dataMat_ = x.visualizations.pca;
	}

	var dataMat = Array();
	n = dataMat_.x1.length;
	for (var i = 0; i < n; i++) 
	{
		dataMat[i] = Array(dataMat_.x1[i], dataMat_.x2[i]);
	}

	var labels;
	var greyedOut = new Array();
	var sixteenS = new Array();
	var outliers = new Array();
	var additionalInfo = '';

	if (highlight16S)
	{
		sixteenS = new Array(dataMat.length);
		for (var i = 0; i < sixteenS.length; ++i) { sixteenS[i] = false; }

		for (var i in x.rnammer.sixteen_s_per_point)
		{
			sixteenS[i] = true;	
		}
	}

	if (selectedLabels === 'fasta')
	{
		labels = Array.apply(null, Array(dataMat.length)).map(Number.prototype.valueOf, -1); // no labels / black color
		additionalInfo = 'Size: ' + (x.fasta_stats.num_basepairs/1000000).toFixed(2) + ' Mbp, GC-content: ' + (x.fasta_stats.gc_content*100).toFixed(2) + ' %';
	} else if(selectedLabels === 'cc')
	{
		labels = labelsPerPoint(x, x.cluster_estimates.cc.assignments[0].assignment.labels);
		outliers = x.cluster_estimates.cc.assignments[0].assignment.outlierClusters;
	} else if(selectedLabels === 'dip')
	{
		if (selectedReduction === 'dataPca')
		{
			if (!selectedNumClusters)
			{
				selectedNumClusters = x.contamination_state === 'contaminated' ? x.cluster_estimates.validity_pca.estimated_k : 1;
				$('#numClusters' + selectedNumClusters).prop('checked', true);
			}			
			labels = labelsPerPoint(x, x.cluster_estimates.validity_pca.assignments[selectedNumClusters-1].assignment.labels);
			outliers = x.cluster_estimates.validity_pca.assignments[selectedNumClusters-1].assignment.outlierClusters;
		} else if (selectedReduction === 'dataSne')
		{
			if (!selectedNumClusters)
			{
				selectedNumClusters = x.contamination_state === 'contaminated' ? x.cluster_estimates.validity_sne.estimated_k : 1;
				$('#numClusters' + selectedNumClusters).prop('checked', true);
			}
			labels = labelsPerPoint(x, x.cluster_estimates.validity_sne.assignments[selectedNumClusters-1].assignment.labels);
			outliers = x.cluster_estimates.validity_sne.assignments[selectedNumClusters-1].assignment.outlierClusters;
		}
	} else if(selectedLabels === 'kraken')
	{
		labels = labelsPerPoint(x, x.kraken.classification);
		for (var i = 0; i < labels.length; i++)
		{
			if (labels[i] === 'unknown')
			{
				greyedOut.push(true);
			} else
			{
				greyedOut.push(false);
			}
		}
		labels = numericLabels(labels);
		additionalInfo = 'Bacterial background: ' + (x.kraken.bacterial_background*100).toFixed(2) + '%';
	}

    var tooltips = x.visualizations.contig_labels.slice();
    for (var i = 0; i < tooltips.length; i++)
    {
    	var contigIdx = tooltips[i];
    	var contig = x.fasta_stats.included_contigs[contigIdx];    
        var len = x.fasta_stats.contig_lengths[contigIdx];
        var gc = x.fasta_stats.contig_gc[contigIdx];
        tooltips[i] = contig + '<br>Size: ' + (len/1000).toFixed(2) + ' kbp<br>GC-content: ' + (gc*100).toFixed(2) + ' %';
    }
	if(selectedLabels === 'kraken')
	{
		tooltips = labelsPerPoint(x, x.kraken.classification);
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
		sixteenS = removeIndexes(sixteenS, rmIdx);
	}

	showData(0, dataMat, labels, tooltips, greyedOut, sixteenS, width, height, padding);
	updateExport(labels, greyedOut);

	$('#additionalInfo').html(additionalInfo);
}

function export16S(idx)
{
	var str = results[selectedFasta].rnammer.sixteen_s_per_point[idx];
	var blob = new Blob([str], {type: "text/plain;charset=utf-8"});
	saveAs(blob, "export_16s.fasta");	
}

function exportClusterFasta(clusterLabel, fasta, selectedLabels, selectedReduction, selectedNumClusters)
{
	var labelsPerCluster;	
	if(selectedLabels === 'cc')
	{
		labelsPerCluster = results[fasta].cluster_estimates.cc.assignments[0].assignment.labels;
	} else if (selectedLabels === 'dip')
	{
		if(selectedReduction === 'dataPca')
		{
			labelsPerCluster = results[fasta].cluster_estimates.validity_pca.assignments[selectedNumClusters-1].assignment.labels;			
		} else if(selectedReduction === 'dataSne')
		{
			labelsPerCluster = results[fasta].cluster_estimates.validity_sne.assignments[selectedNumClusters-1].assignment.labels;			
		}
	} else if (selectedLabels === 'kraken')
	{
		labelsPerCluster = numericLabels(results[fasta].kraken.classification);
	}

	var contigNames = new Array();
	for (var i = 0; i < results[fasta].fasta_stats.included_contigs.length; i++) 
	{
		if (labelsPerCluster[i] == clusterLabel)
		{
			c = results[fasta].fasta_stats.included_contigs[i];
			contigNames.push(c);
		}
	}	

	var str = "";
	for (var c in contigNames)
	{
		str += ">" + contigNames[c] + "\n";
		str += inputcontigs[contigNames[c]] + "\n";
	}

	var blob = new Blob([str], {type: "text/plain;charset=utf-8"});
	saveAs(blob, "export.fasta");	
}

function updateExport(labels, greyedOut)
{
	if (selectedLabels === 'fasta')
	{
		$('#export').hide();
		return;
	}

	// find unknown label
	var unknownLbl = -1;
	for (var i in labels)
	{
		if (greyedOut[i])
		{
			unknownLbl = labels[i];
			break;
		}
	}

	var uniqueLabels = arrayUnique(labels).sort();

	$('#exportColors').html('');
	for (var i in uniqueLabels)
	{
		var lbl = uniqueLabels[i];
		var color = (lbl == unknownLbl) ? greyColor : colors[lbl % colors.length];
		var dataParams = "{ 'clusterLabel': '" + lbl + "', 'fasta': '" + selectedFasta + "', 'selectedLabels': '" + selectedLabels + "', 'selectedReduction': '" + selectedReduction + "', 'selectedNumClusters': '" + selectedNumClusters + "' }";
		var exportLink = $('<a data-params="' + dataParams + '" style="color: ' + color + ';" href="#">&#x25CF;</a>').click(function() 
		{
			var me = $(this), data = me.data('params');
			data = data.replace(/\'/g, '"');
			var obj = jQuery.parseJSON(data);
    		exportClusterFasta(obj.clusterLabel, obj.fasta, obj.selectedLabels, obj.selectedReduction, obj.selectedNumClusters);
		});
		$('#exportColors').append('<span></span>').append(exportLink);		
	}

	$('#export').show();
}
