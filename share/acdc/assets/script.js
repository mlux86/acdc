var colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#a65628", "#333333", "#ffff33"];
var greyColor = '#b3b3b3';
var highlightColor = '#ff7f00';

var width = 500;
var height = 500;
var padding = 25;

var selectedFasta = '';
var selectedLabels = 'fasta';
var selectedReduction = 'dataSne';
var selectedNumClusters = NaN;
var showOutliers = true;
var highlight16S = true;

var results = Array();

$(document).ready(function() 
{
	var id = 1;
	for(var i in clusterinfo)
	{
		results[i] = jsyaml.load(clusterinfo[i]);
		results[i].id = id;
		id = id + 1;
	}

	buildConfidenceTable();

	$(document).keydown(function(evt) {
		var elem = $('td.active');

		var nextActive;
		switch(evt.which)
		{
			case 37: 
				nextActive = elem.prev();
				evt.preventDefault();
				break;
			case 38: 
				nextActive = elem.closest('tr').prev().children().eq(elem.index());
				evt.preventDefault();		
				break;
			case 39: 
				nextActive = elem.next();
				evt.preventDefault();
				break;
			case 40: 
				nextActive = elem.closest('tr').next().children().eq(elem.index());
				evt.preventDefault();
				break;												
			default: break;
		}

		if (typeof nextActive != 'undefined' && nextActive.hasClass('selectable'))
		{
			elem.removeClass('active');
			nextActive.trigger('click');
		}
		
	});

	$('.dataConf').click(function() {
		$('.numClusters').prop("disabled", true);
		selectedFasta = $(this).text();
		selectedLabels = 'fasta';
		showVisualization();
		setActiveExclusively($(this));
	});

	$('.ccConf').click(function() {
		$('.numClusters').prop("disabled", true);
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'cc';
		showVisualization();
		setActiveExclusively($(this));
	});

	$('.dipConf').click(function() {
		$('.numClusters').prop("disabled", false);
		selectedNumClusters = NaN;
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'dip';
		showVisualization();
		setActiveExclusively($(this));
	});

	$('.kraken').click(function() {
		$('.numClusters').prop("disabled", true);
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'kraken';
		showVisualization();
		setActiveExclusively($(this));
	});

	$('.radioReduction').click(function() {
		selectedReduction = $(this).val();
		showVisualization();
	});

	$('.numClusters').click(function() {
		selectedNumClusters = $(this).val();
		showVisualization();
	});

	$('#showOutliers').click(function() {
		showOutliers = $(this).prop('checked');
		showVisualization();
	});

	$('#highlight16S').click(function() {
		highlight16S = $(this).prop('checked');
		showVisualization();
	});

	selectedFasta = $('#confidences').find('td.dataConf:first').text();
	setActiveExclusively($('#confidences').find('td.dataConf:first'));

	if (!selectedFasta)
	{
		$('#confidences').hide();
		$('#scatterContainer').html('<b>No results to show, yet. Reload page to refresh...<b>');
		return;
	}

	$('.numClusters').prop("disabled", true);
	showVisualization();

    $(window).scroll(stickyScatter);
    stickyScatter();

});
