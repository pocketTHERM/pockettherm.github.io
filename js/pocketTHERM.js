// pop-up mouseover
function updateTooltip(input, index) {
  const tooltips = document.getElementsByClassName('tooltiptext');
    tooltips[index].textContent = input.value;
}

function updateFixed(input, index){
  const tooltips = document.getElementsByClassName('comment');
    tooltips[index].textContent = input.value
}

function validateInput(){
  var input = document.getElementsByClassName('myTextbox');
  var errorLabel = document.getElementsByClassName('errorLabel');
  if (input.value === ''){
    errorLabel.innerHTML = "Please enter the value";
    input.classList.add('error');
  } else {
    errorLabel.innerHTML = '';
    input.classList.remove('error');
  }
}

