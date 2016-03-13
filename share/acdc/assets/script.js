var colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"];

var width = 500;
var height = 500;
var padding = 25;

var selectedFasta = '';
var selectedLabels = 'fasta';
var selectedReduction = 'dataSne';
var selectedData = 'oneshot';

$(document).ready(function() 
{

	buildConfidenceTable(results);

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
		selectedFasta = $(this).text();
		selectedLabels = 'fasta';
		updateBootStrapsSelect();
		showVisualization();
		setActiveExclusively($(this));
	});

	$('.pcaConf').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'pca';
		updateBootStrapsSelect();
		showVisualization();
		setActiveExclusively($(this));
	});
	
	$('.sneConf').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'sne';
		updateBootStrapsSelect();
		showVisualization();
		setActiveExclusively($(this));
	});

	$('.kraken').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'kraken';
		updateBootStrapsSelect();
		showVisualization();
		setActiveExclusively($(this));
	});

	$('.radioReduction').click(function() {
		selectedReduction = $(this).val();
		showVisualization();
	});

	$('#bootstraps').change(function() {
		selectedData = $(this).val();
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

	updateBootStrapsSelect();	
	showVisualization();

    $(window).scroll(stickyScatter);
    stickyScatter();

});
