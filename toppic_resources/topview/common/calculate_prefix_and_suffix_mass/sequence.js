class Sequence {
    prsm_data;

    constructor (prsm_data) {
        this.prsm_data = prsm_data;
    }

    getSequence() {
        let sequence = [];
		let firstposition = this.prsm_data.prsm.annotated_protein.annotation.first_residue_position;
		let lastposition = this.prsm_data.prsm.annotated_protein.annotation.last_residue_position;
		this.prsm_data.prsm.annotated_protein.annotation.residue.forEach(function(eachrow){
			if(parseInt(eachrow.position) >= parseInt(firstposition)&&
				parseInt(eachrow.position) <= parseInt(lastposition))
			{
				sequence = sequence+eachrow.acid;
			}
		})
	   return sequence;
    }
}